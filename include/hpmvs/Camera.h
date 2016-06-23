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

#ifndef HPMVS_CAMERA_H_
#define HPMVS_CAMERA_H_

#include <iostream>
#include <hpmvs/NVMReader.h>

namespace mo3d {

class Camera {
public:
	Camera();
	virtual ~Camera();

	void init(const mo3d::NVM_Camera* cam, int width, int height, const int maxLevel = 1);

	/**
	 * Project a 3d world point in homogeneous coordinates into the camera space and
	 * also checks if the point is behind the camera. In addition the resulting point
	 * is clipped by min and max integer values and converted to float.
	 *
	 * @param coord
	 * @param level
	 * @return
	 */
	inline Eigen::Vector3f project(const Eigen::Vector4f& coord, const int level) const {
		if (level >= projection_.size()){
			std::cerr << "Illegal level access in Camera";
			exit(1);
		}
		Eigen::Vector3f result = projection_[level] * coord;

		if (result(2) <= 0.0f) {
			result << -0xffff, -0xffff, -1.0f;
		} else {
			result /= result(2);
			result(0) = std::max((float) (INT_MIN + 3.0f),
					std::min((float) (INT_MAX - 3.0f), result(0)));
			result(1) = std::max((float) (INT_MIN + 3.0f),
					std::min((float) (INT_MAX - 3.0f), result(1)));
		}
		return result;
	}

	inline Eigen::Vector3f project3(const Eigen::Vector3f& coord, const int level) const {
		Eigen::Vector4f hcoord(coord.homogeneous());
		return project(hcoord,level);
	}

	/**
	 * Projects a 3D worldpoint into the camera space, without any additional checks
	 *
	 * @param coord
	 * @param level
	 * @return
	 */
	inline Eigen::Vector3f mult(const Eigen::Vector4f& coord, const int level) const {
		return projection_[level] * coord;
	}

	int getLevels() const { return projection_.size(); }
	float getScale(const Eigen::Vector4f& coord, const int level) const;
	float getLevel(const Eigen::Vector4f& coord, const float scale) const;
	int getLeveli(const Eigen::Vector4f& coord, const float scale, const int maxLevel) const;

public:

	std::string name_;

	// projection matrix for every level
	std::vector<Eigen::Matrix<float, 3, 4> > projection_;

	// internal calibration for every level
	std::vector<Eigen::Matrix<float, 3, 3> > kMat_;

	// position of the camera (optical center)
	Eigen::Vector4f center_;

	// optical axis
	Eigen::Vector4f oAxis_;

	// axis of the image
	Eigen::Vector3f xAxis_, yAxis_, zAxis_;

	// average pixel scale
	float ipscale_;
};

} /* namespace mo3d */

#endif /* HPMVS_CAMERA_H_ */
