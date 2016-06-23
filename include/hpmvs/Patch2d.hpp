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

#ifndef HPMVS_PATCH2D_HPP_
#define HPMVS_PATCH2D_HPP_

#include <Eigen/Dense>

namespace mo3d {

template<unsigned int SIZE, typename PixelT = float>
class Patch2d {
public:
	alignas(128) PixelT data[SIZE*SIZE*3];

	unsigned int size() const {
		return SIZE;
	}

	PixelT dot(const Patch2d<SIZE, PixelT>& other) const {
		const int size = SIZE * SIZE * 3;
		PixelT ans = PixelT(0.0);
		for (size_t i = 0; i < size; ++i) {
			ans += data[i] * other.data[i];
		}
		return ans / size;
	}

	void normalize() {
		const int size = SIZE * SIZE * 3;
		const int size3 = SIZE * SIZE;
		Eigen::Matrix<PixelT, 3, 1> ave = Eigen::Matrix<PixelT, 3, 1>::Zero();

		PixelT* texp = data - 1;
		for (int i = 0; i < size3; ++i) {
			ave[0] += *(++texp);
			ave[1] += *(++texp);
			ave[2] += *(++texp);
		}

		ave /= size3;

		PixelT ave2 = PixelT(0);
		texp = data - 1;
		for (int i = 0; i < size3; ++i) {
			const PixelT f0 = ave[0] - *(++texp);
			const PixelT f1 = ave[1] - *(++texp);
			const PixelT f2 = ave[2] - *(++texp);

			ave2 += f0 * f0 + f1 * f1 + f2 * f2;
		}

		ave2 = std::sqrt(ave2 / size);

		if (ave2 == PixelT(0.0))
			ave2 = PixelT(1.0);

		texp = data - 1;
		for (int i = 0; i < size3; ++i) {
			*(++texp) -= ave[0];
			*texp /= ave2;
			*(++texp) -= ave[1];
			*texp /= ave2;
			*(++texp) -= ave[2];
			*texp /= ave2;
		}
	}

};

typedef Patch2d<7, float> PatchTex;

}

#endif /* HPMVS_PATCH2D_HPP_ */
