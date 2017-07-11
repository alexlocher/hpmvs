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

#ifndef HPMVS_PINHOLE_INTRINSICS_H_
#define HPMVS_PINHOLE_INTRINSICS_H_

namespace mo3d {

struct PinholeIntrinsics {
	std::string imgName;
	int width;
	int height;

	double focal;
	double aspectRatio;
	double skew;

	double ppX;
	double ppY;
	double r1;
	double r2;

	PinholeIntrinsics(const char* name = "", int w = 0, int h = 0, double f = 1.0, double a = 1.0, double s = 0.0,
			double cx = 0.0, double cy = 0.0, double rad1 = 0.0, double rad2 = 0.0) :
			imgName(name), width(w), height(h), focal(f), aspectRatio(a), skew(s), ppX(cx), ppY(cy), r1(rad1), r2(rad2) {
	}

	friend std::ostream& operator <<(std::ostream& os, const PinholeIntrinsics& i) {
		os << i.imgName << " " << i.width << " " << i.height;
		os << " " << i.focal << " " << i.aspectRatio << " " << i.skew;
		os << " " << i.ppX << " " << i.ppY << " " << i.r1 << " " << i.r2;
		os << std::endl;
		return os;
	}

	friend std::istream& operator >>(std::istream& is, PinholeIntrinsics& i){
		is >> i.imgName >> i.width >> i.height;
		is >> i.focal >> i.aspectRatio >> i.skew;
		is >> i.ppX >> i.ppY >> i.r1 >> i.r2;
		return is;
	}
};

}

#endif // HPMVS_PINHOLE_INTRINSICS_H_
