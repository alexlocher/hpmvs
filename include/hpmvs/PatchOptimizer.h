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

#ifndef HPMVS_PATCHOPTIMIZER_H
#define HPMVS_PATCHOPTIMIZER_H

#include <vector>

#include <Eigen/Dense>

#include <hpmvs/Patch3d.h>
#include <hpmvs/Patch2d.hpp>

namespace mo3d {

// forward declaration
class Image;
class Camera;
class HpmvsOptions;
class Scene;

class PatchOptimizer {
public:
	PatchOptimizer(const mo3d::HpmvsOptions& options,const mo3d::Scene* scene);

	bool optimize(mo3d::Patch3d& patch);


private:

	// the optimization fn
	bool runOptimization();

	// image manipulation
	bool addImages();

	void setRefImage();

	bool filterImagesByAngle();

	bool assureImageAngles();

	bool sortImages();

	bool filterImagesNCC(const float threshold);

	// optimization
	bool optimizePatch();

	void setOptimizationFields();

	void getAngleWeightedScales(std::vector<int>& indexes, std::vector<float>& scales,
			std::vector<Eigen::Vector4f>& rays);

	void calculatePatchAxis(const int refIndex, const Eigen::Vector4f& c, const Eigen::Vector4f& n,
			const float scale);

	void setINCCs(std::vector<float> & nccs, const std::vector<int>& indexes, const int refIdx,
			const int robust);

	bool sampleTexture(const Eigen::Vector4f& pCenter, const float pScale,
			const Eigen::Vector4f& pxaxis, const Eigen::Vector4f& pyaxis,
			const Eigen::Vector4f& pzaxis, const int camIdx, mo3d::PatchTex& texture) const;

	void setCenterNorm(const double* x);

	void parametersFromCenterNorm(const Eigen::Vector4f& c, const Eigen::Vector4f n,
			const std::vector<double>& lb, const std::vector<double>& ub, double* x);

	double objective_fn();

	static double static_objective_fn(unsigned n, const double *x, double *grad,
			void *my_func_data);

	static inline float robustincc(const float rhs) {
		return rhs / (1 + 3 * rhs);
	}

	static inline float unrobustincc(const float rhs) {
		return rhs / (1 - 3 * rhs);
	}

private:

	// the optimization memory => we copy the stuff to improve speed
	Eigen::Vector4f pCenter_, pNormal_;
	float pScale_;
	Eigen::Vector4f pXaxis_, pYaxis_, pZaxis_;
	std::vector<int> pImages_;

	// reference related stuff:
	Eigen::Vector4f refCenter_;
	Eigen::Vector4f refRay_;

	// axes of the images
	std::vector<Eigen::Vector3f> imgX_;
	std::vector<Eigen::Vector3f> imgY_;
	std::vector<Eigen::Vector3f> imgZ_;

	float depthScale_;
	float angleScale_;

	// textures
	mo3d::PatchTex refTex_;
	mo3d::PatchTex comTex_;

	// the options
	const mo3d::HpmvsOptions* options_p;

	const mo3d::Image* images_p;
	const mo3d::Camera* camera_p;
	const std::vector<std::vector<int>>* covis_p;
	const mo3d::Scene* scene_p;


};

}

#endif // HPMVS_PATCHOPTIMIZER_H
