/**
* This file is part of HPMVS (Hierarchical Prioritized Multiview Stereo).
*
* Copyright (C) 2015-2016 Alex Locher <alocher at ethz dot ch> (ETH Zuerich)
* For more information see <https://github.com/alexlocher/hpmvs>
*
* HPMVS is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* HPMVS is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with HPMVS. If not, see <http://www.gnu.org/licenses/>.
*/

#include <numeric>
#include <set>
#include <iomanip>

// 3rd party
#include <nlopt.hpp>
#include <glog/logging.h>


// hpmvs stuff
#include <hpmvs/HpmvsOptions.h>
#include <hpmvs/PatchOptimizer.h>
#include <hpmvs/Scene.h>

namespace mo3d {


PatchOptimizer::PatchOptimizer(const mo3d::HpmvsOptions& options, const mo3d::Scene* scene) : options_p(&options), scene_p(scene) {
		camera_p = scene->cameras_.data();
	images_p = scene->images_.data();
	covis_p = &scene->covis_;

	VLOG(2) << "initialized";

}


bool PatchOptimizer::runOptimization() {
	// prepare the patches images
	if (!addImages())
		return false;
	if (!filterImagesNCC(options_p->NCC_ALPHA_1))
		return false;
	sortImages();
	if (!assureImageAngles())
		return false;

	// now do the optimization
	if (!optimizePatch())
		return false;

	// re-handle the attached images for further processing:
	if (!addImages())
		return false;
	if (!filterImagesNCC(options_p->NCC_ALPHA_2))
		return false;
	if (!filterImagesByAngle())
		return false;
	if (!assureImageAngles())
		return false;
	setRefImage();
	if (!filterImagesNCC(options_p->NCC_ALPHA_2))
		return false;

	return true;
}

bool PatchOptimizer::optimize(mo3d::Patch3d& patch) {

	// copy relevant stuff to fields
	this->pCenter_ = patch.center_;
	this->pNormal_ = patch.normal_;
	this->pScale_ = patch.scale_3dx_;
	this->pImages_ = patch.images_;

	if (!runOptimization())
		return false;

	// ok, the optimization was fine => save the values back
	patch.center_ = pCenter_;
	patch.normal_ = pNormal_;
	patch.scale_3dx_ = pScale_;
	patch.images_ = pImages_;

	patch.ncc_ = 1.4f; // FIXME

//	patch. = patch.score2(options_p->NCC_ALPHA_2);

	// set the color
	patch.color_ = scene_p->getColor(patch);

	return true;
}

bool PatchOptimizer::assureImageAngles() {
	// make sure that at least one pair of images fulfills min max constraint

	// get the rays
	std::vector<Eigen::Vector4f> rays;
	for (const int img : pImages_)
		rays.emplace_back((camera_p[img].center_ - pCenter_).normalized());

	// check angle of ray pairs
	const int nrImgs = pImages_.size();
	for (int ii = 0; ii < nrImgs - 1; ii++) {
		for (int jj = ii + 1; jj < nrImgs; jj++) {
			const float a = std::acos(rays[ii].dot(rays[jj]));
			if (a < options_p->MAX_ANGLE && a > options_p->MIN_ANGLE)
				return true;
		}
	}
	return false;
}

bool PatchOptimizer::filterImagesByAngle() {
	std::vector<int> newImages;
	for (const int imgId : pImages_) {
		if ((camera_p[imgId].center_ - pCenter_).normalized().dot(pNormal_)
				> std::cos(options_p->MAX_ANGLE)) {
			newImages.push_back(imgId);
		}
	}
	pImages_.swap(newImages);
	return pImages_.size() >= options_p->MIN_IMAGES_PER_PATCH;

}

bool PatchOptimizer::filterImagesNCC(const float threshold) {
	std::vector<float> inccs;
	setINCCs(inccs, pImages_, 0, 0); // idx 0 as reference and not robustified...

	//----------------------------------------------------------------------
	// Constraint images
	std::vector<int> newimages;
	newimages.push_back(pImages_[0]);
	for (int i = 1; i < (int) pImages_.size(); ++i) {
		if (inccs[i] < 1.0f - threshold)
			newimages.push_back(pImages_[i]);
	}
	pImages_.swap(newimages);
	return pImages_.size() >= options_p->MIN_IMAGES_PER_PATCH;
}

void PatchOptimizer::setRefImage() {
	if (pImages_.size() <= 1)
		return;

	// compute average incc of all pairs
	std::vector<float> incc;
	int refindex = -1;
	float refncc = std::numeric_limits<float>::max();
	for (int ii = 0; ii < pImages_.size(); ii++) {
		setINCCs(incc, pImages_, ii, 1);
		const float sum = accumulate(incc.begin(), incc.end(), 0.0f);
		if (sum < refncc) {
			refncc = sum;
			refindex = ii;
		}
	}

	// swap the indexes
	const int refIndex = pImages_[refindex];
	for (int i = 0; i < (int) pImages_.size(); ++i) {
		if (pImages_[i] == refIndex) {
			const int itmp = pImages_[0];
			pImages_[0] = refIndex;
			pImages_[i] = itmp;
			break;
		}
	}
}

bool PatchOptimizer::sortImages() {
	const float threshold = 1.0f - cos(10.0 * M_PI / 180.0);
	std::vector<int> indexes, indexes2;
	std::vector<float> wScales, wScales2;
	std::vector<Eigen::Vector4f> rays, rays2;

	getAngleWeightedScales(indexes, wScales, rays);

	pImages_.clear();
	if (indexes.size() < 2)
		return false;

	wScales[0] = 0.0f; // keep the reference image!

	while (!indexes.empty()) {
		std::vector<float>::iterator ite = min_element(wScales.begin(), wScales.end());
		const int index = ite - wScales.begin();

		pImages_.push_back(indexes[index]);

		// Remove other images within 5 degrees
		indexes2.clear();
		wScales2.clear();
		rays2.clear();
		for (int j = 0; j < (int) rays.size(); ++j) {
			if (j == index)
				continue;

			indexes2.push_back(indexes[j]);
			rays2.push_back(rays[j]);
			const float ftmp = std::min(threshold,
					std::max(threshold / 2.0f, 1.0f - rays[index].dot(rays[j])));

			wScales2.push_back(wScales[j] * (threshold / ftmp));
		}
		indexes2.swap(indexes);
		wScales2.swap(wScales);
		rays2.swap(rays);
	}
	return pImages_.size() >= options_p->MIN_IMAGES_PER_PATCH;
}

bool PatchOptimizer::addImages() {

	if (pImages_.size() <= 0)
		return false;
	const int refImg = pImages_[0];

	// make sure we do not add existing images again:
	std::set<int> existing_imgs(pImages_.begin(), pImages_.end());
	for (const int covisImg : (*covis_p)[refImg]) {
		if (existing_imgs.find(covisImg) != existing_imgs.end())
			continue;

		// check the angle between normal and ray
		if ((camera_p[covisImg].center_ - pCenter_).normalized().dot(pNormal_)
				< std::cos(options_p->MAX_ANGLE))
			continue;

		// make sure we don't respect to small images
		int imgLevel = std::round(camera_p[covisImg].getLevel(pCenter_, pScale_));
		if (imgLevel < 0 || imgLevel >= options_p->MAXLEVEL - 2)
			continue;

		// check the patche's visibility in the image
		Eigen::Vector2f imgC = camera_p[covisImg].project(pCenter_, imgLevel).head(2);
		if (imgC[0] < 0.0f || images_p[covisImg].getWidth(imgLevel) - 1 <= imgC[0] || imgC[1] < 0.0f
				|| images_p[covisImg].getHeight(imgLevel) - 1 <= imgC[1])
			continue;

		// ok, img is fine => add it
		pImages_.push_back(covisImg);
	}

	return pImages_.size() >= options_p->MIN_IMAGES_PER_PATCH;
}

void PatchOptimizer::getAngleWeightedScales(std::vector<int>& indexes, std::vector<float>& scales,
		std::vector<Eigen::Vector4f>& rays) {
	std::vector<int>::const_iterator bimage = pImages_.begin();
	std::vector<int>::const_iterator eimage = pImages_.end();
	if (pImages_.empty())
		return;
	const int refLevel = std::max(0,
			std::min(options_p->MAXLEVEL - 1,
					(int) std::round(camera_p[pImages_[0]].getLevel(pCenter_, pScale_))));

	indexes.clear();
	scales.clear();
	rays.clear();

	for (const int imgIdx : pImages_) {
		Eigen::Vector4f ray((camera_p[imgIdx].center_ - pCenter_).normalized());
		float cosa = ray.dot(pNormal_.normalized());
		if (cosa > 0) {
			indexes.emplace_back(imgIdx);
			rays.emplace_back(ray);
			float scale = camera_p[imgIdx].getScale(pCenter_, refLevel);
			scales.emplace_back(scale / cosa);
		}
	}
}

double PatchOptimizer::objective_fn() {

	// calculate the patche's axes
	calculatePatchAxis(pImages_[0], pCenter_, pNormal_, pScale_);

	// get the refernce texture
	if (!sampleTexture(pCenter_, pScale_, pXaxis_, pYaxis_, pZaxis_, pImages_[0], refTex_))
		return 2.0;

	// calculate the average incc
	double val = 0.0;
	int nImgs = 0;
	const int optImgs = std::min((int) pImages_.size(), options_p->MAX_IMAGES_PER_PATCH);
	for (int ii = 1; ii < pImages_.size(); ii++) {
		if (!sampleTexture(pCenter_, pScale_, pXaxis_, pYaxis_, pZaxis_, pImages_[ii], comTex_))
			continue;

		val += robustincc(1.0 - refTex_.dot(comTex_));
		nImgs++;
	}
	//if (denom < m_one->options_p->MIN_IMAGES_PER_PATCH - 1)
	if (nImgs < options_p->MIN_IMAGES_PER_PATCH - 1)
		return 2.0;
	else
		return val / nImgs;
}

double PatchOptimizer::static_objective_fn(unsigned n, const double *x, double *grad, void *my_func_data) {

	// recover the pointer tho the PatchOptimizer object
	PatchOptimizer* optimizer = static_cast<PatchOptimizer*>(my_func_data);
	optimizer->setCenterNorm(x);

	return optimizer->objective_fn();
}

bool PatchOptimizer::optimizePatch() {
	if (pImages_.size() < options_p->MIN_IMAGES_PER_PATCH)
		return false;

	double min_angle = -23.99999;	//(- M_PI / 2.0) / m_one->m_ascalesT[id];
	double max_angle = 23.99999;	//(M_PI / 2.0) / m_one->m_ascalesT[id];

	std::vector<double> lower_bounds(3);
	lower_bounds[0] = -HUGE_VAL;		// Not bound
	lower_bounds[1] = min_angle;
	lower_bounds[2] = min_angle;
	std::vector<double> upper_bounds(3);
	upper_bounds[0] = HUGE_VAL;		// Not bound
	upper_bounds[1] = max_angle;
	upper_bounds[2] = max_angle;

	bool success = false;
	double lastVal = 0;
	std::vector<double> x(3, 0);

	try {
		// LN_NELDERMEAD: Corresponds to the N-Simplex-Algorithm of GSL, that was used originally here
		// LN_SBPLX
		// LN_COBYLA
		// LN_BOBYQA
		// LN_PRAXIS
		nlopt::opt opt(nlopt::LN_BOBYQA, 3);
		opt.set_min_objective(static_objective_fn, this);
		opt.set_xtol_rel(1.e-7);
//    opt.set_xtol_rel(1.e-5);
		opt.set_maxeval(1000);

		opt.set_lower_bounds(lower_bounds);
		opt.set_upper_bounds(upper_bounds);

		setOptimizationFields();

		// set initial parameters:
		parametersFromCenterNorm(refCenter_, pNormal_, lower_bounds, upper_bounds, x.data());

//    iterations = 0;
		double minf;
		nlopt::result result = opt.optimize(x, minf);
		lastVal = opt.last_optimum_value();

		success = (result == nlopt::SUCCESS || result == nlopt::STOPVAL_REACHED
				|| result == nlopt::FTOL_REACHED || result == nlopt::XTOL_REACHED);
	} catch (std::exception &e) {
		LOG(WARNING) <<  "patch optimization failed: <" << e.what() << ">";
		success = false;
	}

	if (success) {
		setCenterNorm(x.data());
//		patch.m_ncc = 1.0 - unrobustincc(lastVal); -> FIXME provide last value to fields
	} else {
		return false;
	}

	return true;
}

void PatchOptimizer::setOptimizationFields() {
	imgX_.clear();
	imgY_.clear();
	imgZ_.clear();
	for (int img : pImages_) {
		imgX_.push_back(camera_p[img].xAxis_.normalized());
		imgY_.push_back(camera_p[img].yAxis_.normalized());
		imgZ_.push_back(camera_p[img].zAxis_.normalized());
	}

	refCenter_ = pCenter_;
	refRay_ = (refCenter_ - camera_p[pImages_[0]].center_).normalized();

	depthScale_ = 1.0;
	angleScale_ = M_PI / 48.0f;
}

void PatchOptimizer::setCenterNorm(const double* x) {
	// get center from depth
	pCenter_ = refCenter_ + x[0] * refRay_ * depthScale_;

	// normal from the angles
	const float angle1 = x[1] * angleScale_;
	const float angle2 = x[2] * angleScale_;

	const float fx = sin(angle1) * cos(angle2);
	const float fy = sin(angle2);
	const float fz = -cos(angle1) * cos(angle2);

	pNormal_ << imgX_[0] * fx + imgY_[0] * fy + imgZ_[0] * fz, 0.0;
}

void PatchOptimizer::parametersFromCenterNorm(const Eigen::Vector4f& c, const Eigen::Vector4f n,
		const std::vector<double>& lb, const std::vector<double>& ub, double* x) {
	// depth:
	x[0] = (c - refCenter_).dot(refRay_) / depthScale_;

	// angles
	const float fx = imgX_[0].dot(n.head(3));
	const float fy = imgY_[0].dot(n.head(3));
	const float fz = imgZ_[0].dot(n.head(3));

	x[2] = std::asin(fy);
	const float cosb = std::cos(std::max(-1.0, std::min(1.0, x[2])));

	if (cosb == 0.0)
		x[1] = 0.0;
	else {
		const double sina = fx / cosb;
		const double cosa = -fz / cosb;
		x[1] = std::acos(std::min(1.0, std::max(-1.0, cosa)));
		if (sina < 0.0)
			x[1] = -x[1];
	}

	// scale
	x[1] /= angleScale_;
	x[2] /= angleScale_;

	x[0] = std::min(ub[0], std::max(lb[0], x[0]));
	x[1] = std::min(ub[1], std::max(lb[1], x[1]));
	x[2] = std::min(ub[2], std::max(lb[2], x[2]));
}

void PatchOptimizer::setINCCs(std::vector<float> & inccs, const std::vector<int>& indexes, const int refIdx,
		const int robust) {

	inccs.resize(indexes.size());

	calculatePatchAxis(indexes[refIdx], pCenter_, pNormal_, pScale_);

	// get the reference patch texture
	if (!sampleTexture(pCenter_, pScale_, pXaxis_, pYaxis_, pNormal_, indexes[refIdx], refTex_)) {
		fill(inccs.begin(), inccs.end(), 2.0f);
		return;
	}

	// sample all others on the fly...
	for (int ii = 0; ii < indexes.size(); ii++) {
		if (ii == refIdx)
			inccs[ii] = 0.0f;
		else if (!sampleTexture(pCenter_, pScale_, pXaxis_, pYaxis_, pNormal_, indexes[ii],
				comTex_))
			inccs[ii] = 2.0f;
		else if (robust)
			inccs[ii] = robustincc(1.0f - refTex_.dot(comTex_));
		else
			inccs[ii] = 1.0f - refTex_.dot(comTex_);
	}

}

bool PatchOptimizer::sampleTexture(const Eigen::Vector4f& pCenter, const float pScale,
		const Eigen::Vector4f& pxaxis, const Eigen::Vector4f& pyaxis, const Eigen::Vector4f& pzaxis,
		const int camIdx, mo3d::PatchTex& texture) const {

	const mo3d::Image& image = images_p[camIdx];
	const mo3d::Camera& camera = camera_p[camIdx];

	// check the angle
	if ((camera.center_ - pCenter).normalized().dot(pzaxis) < cos(options_p->MAX_ANGLE))
		return false;

	const int imgLevel = camera.getLeveli(pCenter, pScale, options_p->MAXLEVEL - 1);

	// project the point into the image space
	Eigen::Vector2f center = camera.project(pCenter, imgLevel).head(2);
	Eigen::Vector2f dx = camera.project(pCenter + pxaxis, imgLevel).head(2) - center;
	Eigen::Vector2f dy = camera.project(pCenter + pyaxis, imgLevel).head(2) - center;

	// check the boundaries
	float halfSize = texture.size() / 2.0f;
	Eigen::Vector2f tl = center - halfSize * dx - halfSize * dy; // top left
	Eigen::Vector2f tr = center + halfSize * dx - halfSize * dy; // top right
	Eigen::Vector2f bl = center - halfSize * dx + halfSize * dy; // bottom left
	Eigen::Vector2f br = center + halfSize * dx + halfSize * dy; // bottom right
	Eigen::Vector2f min = tl.cwiseMin(tr).cwiseMin(bl).cwiseMin(br);
	Eigen::Vector2f max = tl.cwiseMax(tr).cwiseMax(bl).cwiseMax(br);

	const int m = 3; // safety margin
	if (min(0) < m || min(1) < m || max(0) >= image.getWidth(imgLevel) - m
			|| max(1) >= image.getHeight(imgLevel) - m) {
		return false;
	}

	// sample for real
	float* target = texture.data;
	Eigen::Vector2f l = tl;
	for (int yy = 0; yy < texture.size(); yy++) {
		Eigen::Vector2f c = l;
		l += dy;
		for (int xx = 0; xx < texture.size(); xx++) {
			Eigen::Vector3f color = scene_p->getColor(camIdx, c[0], c[1], imgLevel);
//			Vec3f color = scene_p->getColor_old(camIdx, c[0], c[1], imgLevel);
//			Vec3f color = m_fm.m_pss.getColor(camIdx, c[0], c[1], imgLevel);
			*(target++) = color[0];
			*(target++) = color[1];
			*(target++) = color[2];
			c += dx;
		}
	}

	texture.normalize();

	return true;
}


void PatchOptimizer::calculatePatchAxis(const int refIndex, const Eigen::Vector4f& c,
		const Eigen::Vector4f& n, const float scale) {
	const mo3d::Camera& refCamera = camera_p[refIndex]; //
	Eigen::Vector3f x, y, z;
	z = n.head(3).normalized();
	y = z.cross(refCamera.xAxis_).normalized();
	x = y.cross(z).normalized();

	x *= scale;
	y *= scale;

	y = y * y.normalized().dot(refCamera.yAxis_.normalized());

	pXaxis_ << x, 0;
	pYaxis_ << y, 0;
	pZaxis_ << z, 0;
}

} // NAMESPACE mo3d
