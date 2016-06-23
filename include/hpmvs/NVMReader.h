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

#ifndef HPMVS_NVMREADER_H_
#define HPMVS_NVMREADER_H_

#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Dense>

namespace mo3d {

struct NVM_Measurement {
	int imgIndex;
	int featIndex;
	Eigen::Vector2d xy;
};

struct NVM_Point {
	Eigen::Vector3d xyz;
	Eigen::Vector3d rgb;
	std::vector<NVM_Measurement> measurements;
};

struct NVM_Camera{
	std::string filename;
	double f; // focal length
	Eigen::Vector4d rq; // rotation quaternion <wxyz>
	Eigen::Vector3d c; // camera center
	double r; // radial distortion
};

struct NVM_Model{
	std::vector<NVM_Camera> cameras;
	std::vector<NVM_Point> points;
};



class NVMReader {
public:
	static void readFile(const char* path, std::vector<NVM_Model>& models, bool fixPath = false);
	static void saveNVM(const char* path, const std::vector<NVM_Model>& models);

	static void saveAsPly(const NVM_Model& model, const char* file);

};

std::istream& operator >>(std::istream& istr, NVM_Measurement& rhs);
std::ostream& operator <<(std::ostream& ostr, NVM_Measurement& rhs);

std::istream& operator >>(std::istream& istr, NVM_Point& rhs);
std::ostream& operator <<(std::ostream& ostr, NVM_Point& rhs);

std::istream& operator >>(std::istream& istr, NVM_Camera& rhs);
std::ostream& operator <<(std::ostream& ostr, NVM_Camera& rhs);

std::istream& operator >>(std::istream& istr, NVM_Model& rhs);
std::ostream& operator <<(std::ostream& ostr, NVM_Model& rhs);

} /* namespace mo3d */

#endif /* HPMVS_NVMREADER_H_ */
