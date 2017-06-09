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

#ifndef HPMVS_TRIANGULATION_HPP_
#define HPMVS_TRIANGULATION_HPP_

#include <Eigen/Dense>

namespace mo3d {

template<typename SCALAR = float>
bool TriangulateMidpoint(const std::vector<Eigen::Matrix<SCALAR, 3, 1>>& origins,
		const std::vector<Eigen::Matrix<SCALAR, 3, 1>>& rays, Eigen::Matrix<SCALAR, 4, 1>* triangulated_point) {
	CHECK_NOTNULL(triangulated_point);
	CHECK_GE(origins.size(), 2);
	CHECK_EQ(origins.size(), rays.size());

	Eigen::Matrix<SCALAR, 4, 4> A;
	A.setZero();
	Eigen::Matrix<SCALAR, 4, 1> b;
	b.setZero();
	for (int i = 0; i < origins.size(); i++) {
		const Eigen::Matrix<SCALAR, 4, 1> ray_direction_homog(rays[i].x(), rays[i].y(), rays[i].z(), 0);
		const Eigen::Matrix<SCALAR, 4, 4> A_term = Eigen::Matrix<SCALAR, 4, 4>::Identity()
				- ray_direction_homog * ray_direction_homog.transpose();
		A += A_term;
		b += A_term * origins[i].homogeneous();
	}

	Eigen::ColPivHouseholderQR<Eigen::Matrix<SCALAR, 4, 4>> qr(A);
	if (qr.info() != Eigen::Success) {
		return false;
	}
	*triangulated_point = qr.solve(b);
	return qr.info() == Eigen::Success;
}


}


#endif /* HPMVS_TRIANGULATION_HPP_ */
