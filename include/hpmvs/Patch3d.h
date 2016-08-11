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

#ifndef HPMVS_PATCH3D_H_
#define HPMVS_PATCH3D_H_

#include <Eigen/Dense>
#include <atomic>
#include <memory>
#include <vector>

namespace mo3d {

class Scene;

class Patch3d {
private:
	static std::atomic_uint_least32_t runningIdx;

		// unique id
		size_t uid_;

		inline uint32_t generateUid() {
			return runningIdx.fetch_add(1);
		}


public:
	Patch3d();
	Patch3d(const Patch3d& other);
	virtual ~Patch3d();
	Patch3d& operator=(const Patch3d &cSource);




	// fields
	Eigen::Vector4f center_;
	Eigen::Vector4f normal_;

	// attached images
	std::vector<int> images_;

	// optimization
	float scale_3dx_;
	float dscale_;

	float ncc_;

	// processor related stuff
	int priorityReduction_;
	bool expanded_;
	float flatness_;

	// color
	Eigen::Vector3f color_;

	// access functions
	inline float x() const {return center_[0];}
	inline float y() const {return center_[1];}
	inline float z() const {return center_[2];}
	inline size_t uid() const {return uid_;}
};

// patch pointer
typedef std::shared_ptr<Patch3d> Ppatch3d;


} /* namespace mo3d */

#endif /* HPMVS_PATCH3D_H_ */
