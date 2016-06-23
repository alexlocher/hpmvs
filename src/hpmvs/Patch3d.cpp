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

#include <hpmvs/Patch3d.h>

namespace mo3d {

// static stuff
std::atomic_uint_least32_t Patch3d::runningIdx(0);

Patch3d::Patch3d() {
	uid_ = generateUid();
	center_ = Eigen::Vector4f::Zero();
	normal_ = Eigen::Vector4f::Zero();
	color_ = Eigen::Vector3f::Zero();

	images_.clear();

	scale_3dx_ = 0.0;
	dscale_ = 0.0;

	ncc_ = 0.0;

	priorityReduction_ = 0;
	expanded_ = false;
	flatness_ = 0.0;
}


Patch3d::Patch3d(const Patch3d& other) {
	uid_ = generateUid();
	*this = other;
}


Patch3d::~Patch3d() {
}



// assignement operator (copy all but uid)
Patch3d& Patch3d::operator=(const Patch3d &cSource) {
	center_ = cSource.center_;
	normal_ = cSource.normal_;
	color_ = cSource.color_;
	images_ = cSource.images_;
	scale_3dx_ = cSource.scale_3dx_;
	dscale_ = cSource.dscale_;
	ncc_ = cSource.ncc_;
	priorityReduction_ = cSource.priorityReduction_;
	expanded_ = cSource.expanded_;
	flatness_ = cSource.flatness_;
	return *this;
}


} /* namespace mo3d */
