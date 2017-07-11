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

#include <hpmvs/Camera.h>
#include <glog/logging.h>

#include <hpmvs/PinholeIntrinsics.h>

namespace mo3d {

Camera::Camera() {
	// TODO Auto-generated constructor stub

}

Camera::~Camera() {
	// TODO Auto-generated destructor stub
}

void Camera::init(const mo3d::NVM_Camera* cam, int width, int height, const int maxLevel){

	projection_.resize(maxLevel+1);
	kMat_.resize(maxLevel+1);

	// calculate projection matrix for level 0
	kMat_[0] <<  cam->f , 0 , width/2.0 , 0 , cam->f , height/2.0 , 0, 0, 1;


	Eigen::Quaterniond q;
	q.w() = cam->rq[0];
	q.x() = cam->rq[1];
	q.y() = cam->rq[2];
	q.z() = cam->rq[3];

	projection_[0].col(3) = -q.matrix().cast<float>() * cam->c.cast<float>() ;
	projection_[0].topLeftCorner(3,3) = q.matrix().cast<float>();
	pMat_ = projection_[0];

	projection_[0] = kMat_[0] * projection_[0];

	// now propagate to lower levels:
	for (int l = 1; l < maxLevel+1; l++){
		projection_[l] = projection_[l-1];
		projection_[l].row(0) /= 2.0;
		projection_[l].row(1) /= 2.0;

		kMat_[l] = kMat_[l-1];
		kMat_[l].row(0) /= 2.0;
		kMat_[l].row(1) /= 2.0;
	}

	// optical center
	center_.head(3) = cam->c.cast<float>();
	center_(3) = 1;

	// axis of the camera
	oAxis_ =  projection_[0].row(2);
	oAxis_ /= projection_[0].row(2).head(3).norm();

	zAxis_ = oAxis_.head(3);
	xAxis_ = projection_[0].row(0).head(3);
	yAxis_ = zAxis_.cross(xAxis_).normalized();
	xAxis_ = yAxis_.cross(zAxis_).normalized();

	// get the average pixelscale
	ipscale_ = (projection_[0].row(0).head(3).norm() + projection_[0].row(1).head(3).norm())/2.0;


	r1_ = cam->r;
	r2_ = 0;

}


void Camera::init(const mo3d::NVM_Camera* cam, const mo3d::PinholeIntrinsics* intrinsics, const int maxLevel){
	CHECK(cam != nullptr && intrinsics != nullptr) << "nullptr as input";

	LOG(INFO)<< "camera >" << cam->filename << "< has focal = " << intrinsics->focal << ", distortion r2 = " << intrinsics->r2;
	projection_.resize(maxLevel + 1);
	kMat_.resize(maxLevel + 1);

	// calculate projection matrix for level 0
	kMat_[0] << intrinsics->focal, intrinsics->skew, intrinsics->ppX, 0, intrinsics->focal * intrinsics->aspectRatio, intrinsics->ppY, 0, 0, 1;

	Eigen::Quaterniond q;
	q.w() = cam->rq[0];
	q.x() = cam->rq[1];
	q.y() = cam->rq[2];
	q.z() = cam->rq[3];

	projection_[0].col(3) = -q.matrix().cast<float>() * cam->c.cast<float>();
	projection_[0].topLeftCorner(3, 3) = q.matrix().cast<float>();
	pMat_ = projection_[0];

	projection_[0] = kMat_[0] * projection_[0];

	// now propagate to lower levels:
	for (int l = 1; l < maxLevel + 1; l++) {
		projection_[l] = projection_[l - 1];
		projection_[l].row(0) /= 2.0;
		projection_[l].row(1) /= 2.0;

		kMat_[l] = kMat_[l - 1];
		kMat_[l].row(0) /= 2.0;
		kMat_[l].row(1) /= 2.0;
	}

	// optical center
	center_.head(3) = cam->c.cast<float>();
	center_(3) = 1;

	// axis of the camera
	oAxis_ = projection_[0].row(2);
	oAxis_ /= projection_[0].row(2).head(3).norm();

	zAxis_ = oAxis_.head(3);
	xAxis_ = projection_[0].row(0).head(3);
	yAxis_ = zAxis_.cross(xAxis_).normalized();
	xAxis_ = yAxis_.cross(zAxis_).normalized();

	// get the average pixelscale
	ipscale_ = (projection_[0].row(0).head(3).norm() + projection_[0].row(1).head(3).norm()) / 2.0;

	r1_ = intrinsics->r1;
	r2_ = intrinsics->r2;

}

float Camera::getScale( const Eigen::Vector4f& coord, const int level) const {
	const float fz = (coord - center_).norm();
	const float ftmp = (kMat_[0](0,0) + kMat_[0](1,1)); // contains fx + fy
	if (ftmp == 0.0)
		return 1.0;

	return 2.0 * fz * (0x0001 << level) / ftmp;
}

float Camera::getLevel(const Eigen::Vector4f& coord, const float scale) const {
	float fz = ( coord - center_).norm();
	return std::log2(scale*(float)(kMat_[0](0,0)+ kMat_[0](1,1))/(2.0 * fz));
}

int Camera::getLeveli(const Eigen::Vector4f& coord, const float scale, const int maxLevel) const {
	return std::max(0, std::min(maxLevel, (int) std::round(getLevel(coord, scale))));
}




} /* namespace mo3d */
