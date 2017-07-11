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

#ifndef HPMVS_IMAGE_IMAGE_H_
#define HPMVS_IMAGE_IMAGE_H_

#include <string>
#include <vector>
#include <Eigen/Dense>


// CImg
#define cimg_display 0
#if defined(PMVS_HAVE_PNG)		// See CMakeLists.txt
#	define cimg_use_png
#endif
#define cimg_use_jpeg
#if defined(PMVS_HAVE_TIFF)		// See CMakeLists.txt
#	define cimg_use_tiff
#endif
//#define cimg_use_zlib
#if defined(_DEBUG)
#	define cimg_verbosity 2
#else
#	define cimg_verbosity 0
#endif
#define WIN32_LEAN_AND_MEAN
#include <CImg.h>



namespace mo3d {

// forward declaration
class NVM_Camera;

class Image {
public:
	Image();
	virtual ~Image();

	void init(const mo3d::NVM_Camera* cam, const int maxLevel = 1);
	bool load(const bool ignoreDistortion = false);

	// confusing because interleaved
	inline int getWidth(int level = 0) const  {return images_[level].height();}
	inline int getHeight(int level = 0) const {return images_[level].depth();}

	inline cimg_library::CImg<unsigned char> getImage(const int level) const {return images_[level].get_permute_axes("yzcx");}

	//inline Eigen::Vector3f getColor(const float fx, const float fy, const int level) const;
	//inline Eigen::Vector3f getColor(const int ix, const int iy, const int level) const;


	inline Eigen::Vector3f getColor(const float x, const float y, const int level) const;

private:
	// the actual image pyramid
	std::vector<cimg_library::CImg<unsigned char> > images_; // storage of the actual images

	// fields needed for loading the image
	std::string path_;
	float f_, k1_;
	int maxLevel_;


	bool undistort(); // undistort Level 0 image

};


Eigen::Vector3f Image::getColor(const float x, const float y, const int level) const{

	// stored interleaved!

  const int lx = static_cast<int>(x);
  const int ly = static_cast<int>(y);
  const int index = 3 * (ly * getWidth(level) + lx);

  const float dx1 = x - lx;  const float dx0 = 1.0f - dx1;
  const float dy1 = y - ly;  const float dy0 = 1.0f - dy1;

  const float f00 = dx0 * dy0;  const float f01 = dx0 * dy1;
  const float f10 = dx1 * dy0;  const float f11 = dx1 * dy1;
  const int index2 = index + 3 * getWidth(level);

  const unsigned char* ucp0 = &images_[level]._data[index] - 1;
  const unsigned char* ucp1 = &images_[level]._data[index2] - 1;
  float r = 0.0f;  float g = 0.0f;  float b = 0.0f;

  r += *(++ucp0) * f00 + *(++ucp1) * f01;
  g += *(++ucp0) * f00 + *(++ucp1) * f01;
  b += *(++ucp0) * f00 + *(++ucp1) * f01;
  r += *(++ucp0) * f10 + *(++ucp1) * f11;
  g += *(++ucp0) * f10 + *(++ucp1) * f11;
  b += *(++ucp0) * f10 + *(++ucp1) * f11;
  return Eigen::Vector3f(r, g, b);
};


} /* namespace mo3d */

#endif /* HPMVS_IMAGE_H_ */
