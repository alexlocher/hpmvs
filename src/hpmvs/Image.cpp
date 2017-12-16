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

#include <hpmvs/Image.h>
#include <glog/logging.h>

#include <nvmtools/NVMFile.h>

namespace mo3d {

Image::Image() :
		path_(""), f_(1.0), k1_(0.0), maxLevel_(0) {
}

Image::~Image() {
}

void Image::init(const nvmtools::NVM_Camera* cam, const int maxLevel) {
	path_ = cam->filename;
	f_ = cam->f;
	k1_ = cam->r;
	maxLevel_ = std::max(1, maxLevel);
}

bool Image::load(const bool ignoreDistortion) {
	images_.clear();
	images_.resize(maxLevel_ + 1);

	// load image on level 0
	images_[0].load(path_.c_str());
	if (images_[0].width() <= 0 || images_[0].height() <= 0)
		return false;

	// check if we have to undistort
	if (k1_ != 0 && ignoreDistortion == false)
		if (!undistort())
			return false;

	// now generate the pyramid
	for (int ii = 1 ; ii < images_.size(); ii++)
		images_[ii] = images_[ii-1].get_resize_halfXY();

	VLOG(2) <<  "loaded pyramid with depth " << images_.size() << " for " << path_;

	// change to interleaved format
	for (int ii = 0; ii < images_.size(); ii++)
		images_[ii].permute_axes("cxyz");

	return true;
}

bool Image::undistort() {
	if (images_.size() < 1 || images_[0].width() <= 0 || images_[0].height() <= 0)
		return false;

	// undistort using the camera model of VisualSFM
	// code from https://groups.google.com/forum/#!msg/vsfm/IcbdIVv_Uek/Us32SBUNK9oJ
	// compute the corresponding distorted pixel coordinate

	const int width = images_[0].width();
	const int height = images_[0].height();

	cimg_library::CImg<float> undistorted(width, height, 1, 3);

	// loop through the image
	for (int ix = 0; ix < width; ix++) {
		for (int iy = 0; iy < height; iy++) {
			float y = (float) (iy - height / 2.0);
			float x = (float) (ix - width / 2.0);

			x /= f_;
			y /= f_;

			if (y == 0)
				y = 1e-3;

			// undistorted coords:
			float mx, my;

			if (k1_ == 0) {
				mx = x;
				my = y;
			} else {
				const double t2 = y * y;
				const double t3 = t2 * t2 * t2;
				const double t4 = x * x;
				const double t7 = k1_ * (t2 + t4);
				if (k1_ > 0) {
					const double t8 = 1.0 / t7;
					const double t10 = t3 / (t7 * t7);
					const double t14 = sqrt(t10 * (0.25 + t8 / 27.0));
					const double t15 = t2 * t8 * y * 0.5;
					const double t17 = pow(t14 + t15, 1.0 / 3.0);
					const double t18 = t17 - t2 * t8 / (t17 * 3.0);
					mx = t18 * x / y;
					my = t18;
				} else {
					const double t9 = t3 / (t7 * t7 * 4.0);
					const double t11 = t3 / (t7 * t7 * t7 * 27.0);
					const std::complex<double> t12 = t9 + t11;
					const std::complex<double> t13 = sqrt(t12);
					const double t14 = t2 / t7;
					const double t15 = t14 * y * 0.5;
					const std::complex<double> t16 = t13 + t15;
					const std::complex<double> t17 = pow(t16, 1.0 / 3.0);
					const std::complex<double> t18 = (t17 + t14 / (t17 * 3.0))
							* std::complex<double>(0.0, sqrt(3.0));
					const std::complex<double> t19 = -0.5 * (t17 + t18) + t14 / (t17 * 6.0);
					mx = t19.real() * x / y;
					my = t19.real();
				}
			}

			x = mx * (float) f_ + width / 2.0f;
			y = my * (float) f_ + height / 2.0f;

			if (x > 1 && x < width - 1 && y > 1 && y < height - 1) {
				undistorted.atXY(ix, iy, 0, 0) = images_[0]._linear_atXY(x, y, 0, 0);
				undistorted.atXY(ix, iy, 0, 1) = images_[0]._linear_atXY(x, y, 0, 1);
				undistorted.atXY(ix, iy, 0, 2) = images_[0]._linear_atXY(x, y, 0, 2);

			}
		}
	}

	// done, swap the data
	images_[0].clear();
	images_[0] = undistorted;

	LOG(INFO) <<  "image undistorted";

	return true;
}

} /* namespace mo3d */
