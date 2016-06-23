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

#include <hpmvs/PatchOptimizer.h>
#include <hpmvs/Scene.h>

#include <chrono>

#include <glog/logging.h>
#include <stlplus3/file_system.hpp>

namespace mo3d {

// static fields
float Scene::MAX_DEPTH = 1000.0f;

Scene::Scene() {
	t_offtime = 0;
}

Scene::~Scene() {
}

bool Scene::addCameras(const NVM_Model& model, const HpmvsOptions& options) {

	// iterate through the cameras and create corresponding objects
	for (const NVM_Camera& cam : model.cameras) {
		int newCamId = cameras_.size();

		// create new objects
		cameras_.emplace_back();
		images_.emplace_back();
		m_depth_lock.emplace_back();
		m_depths.emplace_back();
		dict_[cam.filename] = newCamId;
	}

	// load the images
	std::atomic<int> imgsLoaded(0);
	const int totalImgs = model.cameras.size();

#pragma omp parallel for
	for (int ii = 0; ii < model.cameras.size(); ii++) {

		// get the index/id
		const int camId = dict_[model.cameras[ii].filename];
		CHECK (camId < cameras_.size()) << "we have a problem";

		// populate them
		images_[camId].init(&model.cameras[ii], options.MAXLEVEL);
		images_[camId].load();

		cameras_[camId].init(&model.cameras[ii], images_[camId].getWidth(),
				images_[camId].getHeight(), options.MAXLEVEL);

		// corresponding depth image
		for (int level = 0; level < cameras_[camId].getLevels(); level++) {
			int rows = images_[camId].getHeight(level) / DEPTH_SUBSAMPLE;
			int cols = images_[camId].getWidth(level) / DEPTH_SUBSAMPLE;
			m_depths[camId].emplace_back(new Eigen::MatrixXf(rows, cols));
			*m_depths[camId].back() = Eigen::MatrixXf::Ones(rows, cols) * MAX_DEPTH;
			m_depth_lock[camId].emplace_back(new std::mutex);
		}

		imgsLoaded++;
		LOG(INFO)<< "Loaded " << imgsLoaded.load() << " / " << totalImgs << " images";
	}

	return true;
}

bool Scene::initPatches(const NVM_Model& model, const HpmvsOptions& options) {
	int cSize = 2;

	// create optimizer instances
	std::vector<PatchOptimizer> optimizers;
	for (int ii = 0; ii < omp_get_max_threads(); ii++) {
		optimizers.emplace_back(options, this);
	}

	const auto start = std::chrono::steady_clock::now();
	t_start = std::chrono::steady_clock::now();
	t_offtime = 0;

	std::vector<Ppatch3d> initPatches;

	// loop through the points in the nvm_model
	const size_t nPts = model.points.size();
#pragma omp parallel for
	for (int ii = 0; ii < nPts; ii++) {
		const NVM_Point& pt = model.points[ii];
		Ppatch3d ppatch(new Patch3d);
		ppatch->center_.head(3) = pt.xyz.cast<float>();
		ppatch->center_[3] = 1.0;

		// only accept point with enough images
		if (pt.measurements.size() < options.MIN_IMAGES_PER_PATCH)
			continue;

		Eigen::Vector4f commonCenter = Eigen::Vector4f::Zero();

		// add all measurements to the patch
		for (const NVM_Measurement& m : pt.measurements) {
//			int idx = dict_[m.imgIndex];
			int idx = m.imgIndex;
			if (idx < 0)
				continue;

			// check visibility of that patch:
			Eigen::Vector3f p = cameras_[idx].project(ppatch->center_, options.START_LEVEL);

			const int margin = cSize;
			if (p[0] < margin || p[1] < margin
					|| p[0] >= images_[idx].getWidth(options.START_LEVEL) - margin
					|| p[1] >= images_[idx].getHeight(options.START_LEVEL) - margin) {
				continue;
			}

			ppatch->images_.push_back(idx);
			commonCenter += cameras_[idx].center_;
		}

		if (ppatch->images_.size() < 2)
			continue; // not enough images attached

		commonCenter = commonCenter / (float) ppatch->images_.size();
		ppatch->normal_ = commonCenter - ppatch->center_;
		ppatch->normal_ = cameras_[ppatch->images_[0]].center_ - ppatch->center_; // just first
		ppatch->normal_.normalize();
		ppatch->normal_[3] = 0.0;
		ppatch->scale_3dx_ = cameras_[ppatch->images_[0]].getScale(ppatch->center_,
				options.START_LEVEL);

		// optimization
		const int threadId = omp_get_thread_num();
		if (!optimizers[threadId].optimize(*ppatch))
			continue;

		// make sure we did not moved to far...
		if ((ppatch->center_.head(3) - pt.xyz.cast<float>()).norm() > ppatch->scale_3dx_ * 2)
			continue;

#pragma omp critical
		{
			initPatches.push_back(ppatch);
		}
	}

	const auto end = std::chrono::steady_clock::now();
	LOG(INFO)<< "Initialized from NVM within " << std::chrono::duration<double>(end - start).count()<< " seconds";
	LOG(INFO)<< "created " << initPatches.size()<< " initPatches out of " << model.points.size()<< " points";

	// now fill the created patches into the tree
	// -----------------------------------------------------------------
	Eigen::Vector3f min, max, dist;
	getBoundingBox(initPatches, min, max);
	dist = max - min;
	float width = std::max(dist[0], std::max(dist[1], dist[2]));

	Branch<Ppatch3d>* oldRoot_p = patchTree_.swapRoot(
			new Branch<Ppatch3d>((min + max) / 2.0, width));
	delete oldRoot_p;

	for (auto& p : initPatches) {
		p->scale_3dx_ = std::max(p->scale_3dx_, width / (1 << 10));
		patchTree_.add(p, p->scale_3dx_);
		setDepths(*p, false);
	}

	patchTree_.cellHistogram();

	// debug
	// visualizeDepths("/tmp");

	return true;

}

bool Scene::extractCoVisiblilty(const NVM_Model& model, const HpmvsOptions& options) {

	const char* outVisDataFile = "/tmp/testvis.dat";

	int nCams = model.cameras.size();
	if (nCams != cameras_.size())
		return false;

	// initialize empty nCams x nCams matrix
	Eigen::MatrixXi vis = Eigen::MatrixXi::Zero(nCams, nCams);

	// fill visibility information
	for (const mo3d::NVM_Point& p : model.points) {
		std::vector<int> visIds;
		for (const NVM_Measurement& m : p.measurements) {
			int camId = m.imgIndex;
			visIds.emplace_back(camId);
		}

		for (int ii = 0; ii < visIds.size(); ii++) {
			for (int jj = 0; jj < visIds.size(); jj++) {
				if (ii != jj)
					vis(ii, jj)++;}
			}
		}

		// now consider cameras to have co-visibility if they share at least 50 points
	covis_.clear();
	covis_.resize(nCams);
	for (int ii = 0; ii < vis.rows(); ii++) {
		for (int jj = 0; jj < vis.cols(); jj++) {
			if (vis(ii, jj) >= 50) {
				covis_[ii].emplace_back(jj);
			}
		}
	}

	// for debug purposes we can also store the resulting visdata file
	if (outVisDataFile != 0) {
		std::ofstream of(outVisDataFile);

		of << "VISDATA" << std::endl; // header
		of << nCams << std::endl; // number of cameras
		for (int ii = 0; ii < covis_.size(); ii++) {
			of << ii << " " << covis_[ii].size();
			for (int jj : covis_[ii])
				of << " " << jj;

			of << std::endl;
		}

		of.close();
	}

	LOG(INFO)<< "extracted visibility information from nvm file";

	return true;
}

Eigen::Vector3f Scene::getColor(const Patch3d& patch) const {

	std::vector<Eigen::Vector3f> colors;

	// loop through the visible images
	for (int imgIdx : patch.images_) {
		const Camera& camera = cameras_[imgIdx];
		const int imgLevel = camera.getLeveli(patch.center_, patch.scale_3dx_,
				camera.getLevels() - 1);
		Eigen::Vector2f c = camera.project(patch.center_, imgLevel).head(2);
//		const cimg_library::CImg<float>& img = images_[imgIdx].getImage(imgLevel);

		colors.push_back(images_[imgIdx].getColor(c(0), c(1), imgLevel));
//		colors.back()[0] = img._linear_atXY(c(0), c(1), 0);
//		colors.back()[1] = img._linear_atXY(c(0), c(1), 1);
//		colors.back()[2] = img._linear_atXY(c(0), c(1), 2);
	}

	// now get the median
	std::sort(colors.begin(), colors.end(),
			[] (const Eigen::Vector3f& a, const Eigen::Vector3f& b) {return a.norm() < b.norm();});

	// handle white pixels (e.g. sky)
	if ((colors[colors.size() / 2]).norm() > 250.0)
		return colors.front();
	else
		return colors[colors.size() / 2];
}

/**
 * gets the number of images, having images attach with higher resolution (not yet L0) available
 *
 * @param patch
 * @return
 */
int Scene::getLevelSupport(const Patch3d& patch) {
	int nrImgs = 0;

	for (int imgIdx : patch.images_) {
		if (std::round(cameras_[imgIdx].getLevel(patch.center_, patch.scale_3dx_)) > 0)
			nrImgs++;
	}

	return nrImgs;
}

/**
 * Sets the depth of the associated visible images within the depth into the depth memory
 *
 * @param patch
 */
void Scene::setDepths(const Patch3d& patch, bool subtract) {
	const int factor = (subtract) ? -1 : 1;
	for (int idx : patch.images_) {
		int level = cameras_[idx].getLeveli(patch.center_, patch.scale_3dx_,
				cameras_[idx].getLevels() - 1);

		const Eigen::Vector3f imgC = cameras_[idx].mult(patch.center_, level);

		int x = (int) (imgC[0] / imgC[2] + 0.5) / DEPTH_SUBSAMPLE;
		int y = (int) (imgC[1] / imgC[2] + 0.5) / DEPTH_SUBSAMPLE;
		float d = imgC[2];

		CHECK (d >= 0) << "negative depth detected, you have a problem";

		if (x < 0 || x >= m_depths[idx][level]->cols() || y < 0
				|| y >= m_depths[idx][level]->rows())
			continue;

		std::lock_guard<std::mutex> lock(*m_depth_lock[idx][level]);

		float oldDepth = (*m_depths[idx][level])(y, x);
		if (oldDepth == d && subtract)
			(*m_depths[idx][level])(y, x) = MAX_DEPTH;
		else if (subtract == false && d < oldDepth)
			(*m_depths[idx][level])(y, x) = d;

//		(*m_depths[idx][level])(y, x) += (imgC[2] * factor);
//		(*m_depths_count[idx][level])(y, x) += factor;

	}
}

float Scene::getDetphAtLevel(const int imgIdx, const int xx, const int yy, const int level,
bool lock) {
	// traverse the depth pyramid and collect depth values
	float depth = MAX_DEPTH;

	int x = xx / DEPTH_SUBSAMPLE;
	int y = yy / DEPTH_SUBSAMPLE;

	if (x < 0 || x >= m_depths[imgIdx][level]->cols() || y < 0
			|| y >= m_depths[imgIdx][level]->rows())
		return depth;

	if (lock) {
		std::unique_lock<std::mutex> celllock((*m_depth_lock[imgIdx][level]));
		depth = (*m_depths[imgIdx][level])(y, x);
		celllock.unlock();
	} else {
		depth = (*m_depths[imgIdx][level])(y, x);
	}

	return depth;
}

float Scene::getFullDepth(const int imgIdx, const int xx, const int yy, bool lock) {
	// traverse the depth pyramid and collect depth values
	float depth = MAX_DEPTH;

	int x = xx / DEPTH_SUBSAMPLE;
	int y = yy / DEPTH_SUBSAMPLE;

	const int levels = cameras_[imgIdx].getLevels();

	for (int level = 0; level < levels; level++) {
		if (x < 0 || x >= m_depths[imgIdx][level]->cols() || y < 0
				|| y >= m_depths[imgIdx][level]->rows())
			return depth;
		if (lock) {
			std::unique_lock<std::mutex> celllock((*m_depth_lock[imgIdx][level]));
			depth = std::min(depth, (*m_depths[imgIdx][level])(y, x));
			celllock.unlock();
		} else {
			depth = std::min(depth, (*m_depths[imgIdx][level])(y, x));
		}
		x /= 2;
		y /= 2;
	}

	return depth;

}

void Scene::visualizeDepths(const char* folder) {
	// create a html file with an image table for visualization
	std::ofstream of(stlplus::create_filespec(folder, "overview", "html").c_str(),
			std::ofstream::out);
	of << "<!DOCTYPE html><html><head>";
	of << "<style>table, th, td {border: 1px solid black;border-collapse: collapse;}";
	of << "img { height: auto; width: 100%;}";
	of << "th, td {padding: 5px;text-align: left;}</style>";
	of << "</head><body>";
	of << "<h2>Depth Images:</h2>";
	of << "<table style=\"width:100%\">";

	// get the number levels and images
	int m_maxLevel = cameras_[0].getLevels();
	int m_num = cameras_.size();

	// html heading
	of << "<tr><th>Color</th><th>Combined</th>";
	for (int level = 0; level < m_maxLevel; level++)
		of << "<th>L" << level << "</th>";
	of << "</tr>";

	// save the images into the provided folder with the prefix idx_level.png
	for (int imgidx = 0; imgidx < m_num; imgidx++) {
		of << "<tr>";

		// colored image
		{
			int level = 1;
			std::string filename = stlplus::create_filename(std::to_string(imgidx) + "_col", "jpg");
			images_[imgidx].getImage(level).save(
					stlplus::create_filespec(folder, filename).c_str());
			of << "<td><img src=\"" << filename << "\"/></td>";
		}

		// accumulated image
		{
			cimg_library::CImg<float> img(images_[imgidx].getWidth(0) / DEPTH_SUBSAMPLE,
					images_[imgidx].getHeight(0) / DEPTH_SUBSAMPLE, 1, 1);
			for (int yy = 0; yy < img.height(); yy++) {
				for (int xx = 0; xx < img.width(); xx++) {
					img._atXY(xx, yy) = getFullDepth(imgidx, xx * DEPTH_SUBSAMPLE,
							yy * DEPTH_SUBSAMPLE, false);
					if (img._atXY(xx, yy) == MAX_DEPTH)
						img._atXY(xx, yy) = 0;
				}
			}
			img.normalize(0, 255.0);
			cimg_library::CImg<unsigned char> cmapImg = img.map(
					cimg_library::CImg<unsigned char>::jet_LUT256());
			std::string filename = stlplus::create_filename(std::to_string(imgidx) + "_all", "jpg");
			cmapImg.save(stlplus::create_filespec(folder, filename).c_str());

			of << "<td><img src=\"" << filename << "\"/></td>";

		}

		for (int level = 0; level < m_maxLevel; level++) {
			// create an image
			cimg_library::CImg<float> img(images_[imgidx].getWidth(level) / DEPTH_SUBSAMPLE,
					images_[imgidx].getHeight(level) / DEPTH_SUBSAMPLE, 1, 1);
			for (int yy = 0; yy < img.height(); yy++) {
				for (int xx = 0; xx < img.width(); xx++) {
					img._atXY(xx, yy) = (*m_depths[imgidx][level])(yy, xx);
					if (img._atXY(xx, yy) == MAX_DEPTH)
						img._atXY(xx, yy) = 0;
				}
			}
			img.normalize(0, 255.0);
			cimg_library::CImg<unsigned char> cmapImg = img.map(
					cimg_library::CImg<unsigned char>::jet_LUT256());
			std::string filename = stlplus::create_filename(
					std::to_string(imgidx) + "_" + std::to_string(level), "jpg");
			cmapImg.save(stlplus::create_filespec(folder, filename).c_str());
			of << "<td><img src=\"" << filename << "\"/></td>";
		}
		of << "</tr>";
	}

	of << "</table></body></html>";
	of.close();

}

int Scene::depthTests(const Patch3d& patch, const float margin) {
	int nrVisible = 0;
	for (int imgIdx : patch.images_)
		if (depthTest(patch, imgIdx, margin, true, false))
			++nrVisible;
	return nrVisible;
}

/**
 * Does a depth test with the provided patch in the given image index.
 *
 * @param patch
 * @param image
 * @param margin
 * @return
 */
bool Scene::depthTest(const Patch3d& patch, const int image, const float margin, bool neighbours,
bool viewBlock) {
	// for the visibility test we always work on level 0
	const int level = 0;
	Eigen::Vector3f imgC = cameras_[image].mult(patch.center_, level);
	int ix = (int) (imgC[0] / imgC[2] + 0.5);
	int iy = (int) (imgC[1] / imgC[2] + 0.5);

	if (!neighbours)
		return depthTest(patch, ix, iy, imgC[2], image, margin, viewBlock);

	ix--;
	iy--;

	for (int yy = 0; yy < 3; yy++) {
		for (int xx = 0; xx < 3; xx++) {
			if (!depthTest(patch, ix + xx, iy + yy, imgC[2], image, margin, viewBlock))
				return false;
		}
	}
	return true;

}

bool Scene::depthTest(const Patch3d& patch, const int ix, const int iy, const float depth,
		const int image, const float margin, bool viewBlock) {
	int level = 0;
	if (depth < 0 || ix < 0 || ix >= images_[image].getWidth(level) || iy < 0
			|| iy >= images_[image].getHeight(level))
		return false;

	const float imgDepth = getFullDepth(image, ix, iy, false);
	if (imgDepth >= MAX_DEPTH)
		return (viewBlock) ? false : true; // in case we test for view blocking, this is not a problem in similar depth, we are good to go

	// just do the same as pmvs FIXME
	Eigen::Vector4f ray = (patch.center_ - cameras_[image].center_).normalized();
	const float diff = imgDepth - depth;
	const float factor = std::min(2.0f, 2.0f + ray.dot(patch.normal_));

	if (viewBlock == false) {
		// testing for similar depth // FIXME had csize in here... => I put 2 now
		return abs(diff) < patch.scale_3dx_ * /* m_fm.m_csize  * */margin * factor * 2.0;
	} else {
		// test if the patch lying in front of a confirmed surface (and hence blocking the view)
		return diff > patch.scale_3dx_ * /* m_fm.m_csize * */margin * factor * 2.0;
	}
	return false;
}

int Scene::pixelFreeTests(const Patch3d& patch) {
	int nrFree = 0;
	for (int imgIdx : patch.images_)
		if (pixelFreeTest(patch, imgIdx))
			++nrFree;
	return nrFree;
}

bool Scene::pixelFreeTest(const Patch3d& patch, const int image) {
	int level = (int) std::round(cameras_[image].getLevel(patch.center_, patch.scale_3dx_));
	if (level < 0 || level >= cameras_[image].getLevels())
		return false;

	Eigen::Vector3f imgC = cameras_[image].project(patch.center_, level);
	int ix = (int) (imgC[0] / imgC[2] + 0.5);
	int iy = (int) (imgC[1] / imgC[2] + 0.5);

	if (ix < 0 || ix >= images_[image].getWidth(level) || iy < 0
			|| iy >= images_[image].getHeight(level))
		return false;

	const float imgDepth = getDetphAtLevel(image, ix, iy, level);
	return (imgDepth == MAX_DEPTH);

}

int Scene::viewBlockTest(const Patch3d& patch, const float margin) {
	int nrBlocking = 0;
	// loop through all the image and check if the patch should be visible
	for (int imgIdx = 0; imgIdx < images_.size(); imgIdx++) {

		int level = (int) std::round(cameras_[imgIdx].getLevel(patch.center_, patch.scale_3dx_));
		if (level < 0 || level > cameras_[imgIdx].getLevels() - 1)
			continue;

		// check if the patch is theoretically visible within that image
		Eigen::Vector3f imgC = cameras_[imgIdx].project(patch.center_, level);
		if (imgC[0] < 0 || imgC[0] > images_[imgIdx].getWidth(level) || imgC[1] < 0
				|| imgC[1] > images_[imgIdx].getHeight(level))
			continue;

		/*
		 // angle (30 deg)
		 Eigen::Vector3f n(patch->m_coord[0], patch->m_coord[1], patch->m_coord[2]);
		 Eigen::Vector3f O(m_fm.m_pss.m_photos[imgIdx].m_oaxis[0], m_fm.m_pss.m_photos[imgIdx].m_oaxis[1],
		 m_fm.m_pss.m_photos[imgIdx].m_oaxis[2]);
		 n.normalize();
		 O.normalize();
		 double angle = acos(n.dot(O));

		 if (angle / M_PI * 180 > 45)
		 continue;
		 */

		// and do the test
		if (depthTest(patch, imgIdx, margin, true, true))
			nrBlocking++;
	}

	return nrBlocking;
}

void Scene::saveAsNVM(const char* folder) {
	// create the output folder if needed
	if (!stlplus::folder_exists(folder))
		stlplus::folder_create(folder);

	// save images
	const std::string imgFolder = stlplus::create_filespec(folder, "imgs");
	if (!stlplus::folder_exists(imgFolder))
		stlplus::folder_create(imgFolder);
	std::vector<std::string> imgNames;
	for (int ii = 0; ii < images_.size(); ii++) {
		imgNames.emplace_back(stlplus::create_filespec(imgFolder, std::to_string(ii), ".jpg"));
		images_[ii].getImage(0).save(imgNames.back().c_str());
	}

	// create an NVM Model out of the scene
	std::vector<NVM_Model> models(1);
	for (int ii = 0; ii < images_.size(); ii++) {
		NVM_Camera cam;
		cam.filename = stlplus::create_filespec("imgs", std::to_string(ii), ".jpg");
		cam.f = cameras_[ii].kMat_[0](0, 0);
		cam.c = cameras_[ii].center_.cast<double>().hnormalized();
		cam.r = 0;
		Eigen::Matrix3d m(Eigen::Matrix3d::Identity());
		m.row(0) << cameras_[ii].xAxis_.normalized().cast<double>().transpose();
		m.row(1) << cameras_[ii].yAxis_.normalized().cast<double>().transpose();
		m.row(2) << cameras_[ii].zAxis_.normalized().cast<double>().transpose();
		Eigen::Quaterniond rq(m);
		cam.rq << rq.w(), rq.x(), rq.y(), rq.z();

		// save it
		models.back().cameras.push_back(cam);
	}

	// now loop through the tree
	Leaf_iterator<Ppatch3d> end = patchTree_.end();
	for (Leaf_iterator<Ppatch3d> leaf = patchTree_.begin(); leaf != end; leaf++) {
		std::vector<Ppatch3d>& patches = leaf->data;
		for (const auto& p : patches) {

			// create a point
			NVM_Point np;
			np.rgb = p->color_.cast<double>();
			np.xyz = p->center_.cast<double>().hnormalized();

			// add measurements
			for (int ii = 0; ii < p->images_.size(); ii++) {
				NVM_Measurement m;
				m.featIndex = 0;
				m.imgIndex = p->images_[ii];

				// we have to project:
				m.xy = cameras_[m.imgIndex].project(p->center_, 0).head(2).cast<double>();

				np.measurements.push_back(m);
			}

			models.back().points.push_back(np);
		}
	}

	models.emplace_back();

	// now save the model
	std::string nvmfile = stlplus::create_filespec(folder, "project", "nvm");
	NVMReader::saveNVM(nvmfile.c_str(), models);

}

void Scene::savePMats(const char* file) {
	std::ofstream of(file, std::ofstream::out);
	Eigen::IOFormat fmt(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
	for (int ii = 0; ii < cameras_.size(); ii++) {
		of << cameras_[ii].projection_[0].format(fmt) << std::endl;
	}
	of.close();
}

void Scene::savePoseMats(const char* file) {
	std::ofstream of(file, std::ofstream::out);
	Eigen::IOFormat fmt(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
	for (int ii = 0; ii < cameras_.size(); ii++) {
		Eigen::Matrix3d m(Eigen::Matrix3d::Identity());
		m.row(0) << cameras_[ii].xAxis_.normalized().cast<double>().transpose();
		m.row(1) << cameras_[ii].yAxis_.normalized().cast<double>().transpose();
		m.row(2) << cameras_[ii].zAxis_.normalized().cast<double>().transpose();
		Eigen::Matrix<double, 3, 4> pose;
		pose << m, cameras_[ii].center_.cast<double>().head(3);

		of << pose.format(fmt) << std::endl;
	}
	of.close();

}

} /* namespace mo3d */
