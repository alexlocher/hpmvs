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

#include <hpmvs/NVMReader.h>
#include <glog/logging.h>
#include <stlplus3/portability/portability.hpp>
#include <iomanip>

using namespace std;


namespace mo3d {

std::istream& operator >>(std::istream& istr, NVM_Measurement& rhs) {
	istr >> rhs.imgIndex >> rhs.featIndex >> rhs.xy[0] >> rhs.xy[1];
	return istr;
}
std::ostream& operator <<(std::ostream& ostr, NVM_Measurement& rhs) {
	ostr << " " << rhs.imgIndex << " " << rhs.featIndex << " " << rhs.xy[0] << " " << rhs.xy[1];
	return ostr;
}

std::istream& operator >>(std::istream& istr, NVM_Point& rhs) {
	istr >> rhs.xyz[0] >> rhs.xyz[1] >> rhs.xyz[2] >> rhs.rgb[0] >> rhs.rgb[1] >> rhs.rgb[2];
	int nrMeasurements;
	istr >> nrMeasurements;
	rhs.measurements.resize(nrMeasurements);
	for (size_t ii = 0; ii < nrMeasurements; ii++) {
		istr >> rhs.measurements[ii];
	}
	return istr;
}

std::ostream& operator <<(std::ostream& ostr, NVM_Point& rhs) {
	ostr << rhs.xyz[0] << " " << rhs.xyz[1] << " " << rhs.xyz[2] << " " << ((int) rhs.rgb[0]) << " "
			<< ((int) rhs.rgb[1]) << " " << ((int) rhs.rgb[2]);
	int nrMeasurements = rhs.measurements.size();
	ostr << " " << nrMeasurements;
	for (size_t ii = 0; ii < nrMeasurements; ii++) {
		ostr << rhs.measurements[ii];
	}
	ostr << endl;
	return ostr;
}

std::istream& operator >>(std::istream& istr, NVM_Camera& rhs) {
	istr >> rhs.filename;
	istr >> rhs.f;
	istr >> rhs.rq[0] >> rhs.rq[1] >> rhs.rq[2] >> rhs.rq[3];
	istr >> rhs.c[0] >> rhs.c[1] >> rhs.c[2];
	istr >> rhs.r;
	int check;
	istr >> check;
	assert(check == 0 && "last camera parameter shoud be 0");
	replace(rhs.filename.begin(), rhs.filename.end(), '"', ' ');
	return istr;
}
std::ostream& operator <<(std::ostream& ostr, NVM_Camera& rhs) {
	ostr << rhs.filename << " ";
	ostr << rhs.f << " ";
	ostr << rhs.rq[0] << " " << rhs.rq[1] << " " << rhs.rq[2] << " " << rhs.rq[3] << " ";
	ostr << rhs.c[0] << " " << rhs.c[1] << " " << rhs.c[2] << " ";
	ostr << rhs.r << " ";
	int check = 0;
	ostr << check << endl;
	return ostr;
}

std::istream& operator >>(std::istream& istr, NVM_Model& rhs) {
	int nrCameras = 0, nrPoints = 0;
	istr >> nrCameras;
	rhs.cameras.resize(nrCameras);
	for (size_t ii = 0; ii < nrCameras; ii++)
		istr >> rhs.cameras[ii];

	// handle empty case
	if (nrCameras > 0)
		istr >> nrPoints;
	rhs.points.resize(nrPoints);
	for (size_t ii = 0; ii < nrPoints; ii++)
		istr >> rhs.points[ii];
	return istr;
}
std::ostream& operator <<(std::ostream& ostr, NVM_Model& rhs) {
	int nrCameras = rhs.cameras.size(), nrPoints = rhs.points.size();
	ostr << endl << nrCameras << endl;
	for (size_t ii = 0; ii < nrCameras; ii++)
		ostr << rhs.cameras[ii];

	// handle empty case
	if (nrCameras > 0)
		ostr << endl << nrPoints << endl;
	for (size_t ii = 0; ii < nrPoints; ii++)
		ostr << rhs.points[ii];
	return ostr;
}

void NVMReader::readFile(const char* path, std::vector<NVM_Model>& models, bool fixPath) {
	models.clear();
	std::ifstream infile(path);
	if (!infile.good()) {
		LOG(WARNING)<< "cannot read from <" << path << ">";
		return;
	}

	// extract the folder of this nvm file
	string nvmfolder(stlplus::folder_part(path));

	// check header
	string header;
	infile >> header;
	if (strcasecmp("NVM_V3", header.c_str()) != 0) {
		LOG(WARNING)<< "<" << path << "> is no valid nvm file [TAG = " << header << "]";
		return;
	}

	do {
		models.emplace_back();
		infile >> models.back();

		if (fixPath) {
			for (size_t ii = 0; ii < models.back().cameras.size(); ii++) {
				string name = models.back().cameras[ii].filename;
				if (stlplus::is_relative_path(name))
					models.back().cameras[ii].filename = stlplus::create_filespec(nvmfolder, name);
			}

		}

	} while (infile.good() && models.back().cameras.size() > 0);
	if (models.size() > 0)
		models.pop_back(); // remove empty model at the back

	LOG(INFO)<< "read " << models.size() << " models from <" << path << ">";
	for (size_t ii = 0; ii < models.size(); ii++)
		LOG(INFO)<< "Model " << ii << ": " << models[ii].cameras.size()<< " cameras and " << models[ii].points.size()<< " points";

	}

void NVMReader::saveNVM(const char* path, const std::vector<NVM_Model>& models) {

	std::ofstream outfile(path);
	if (!outfile.good()) {
		LOG(WARNING)<< "cannot write to <" << path << ">";
		return;
	}

	// more precision
	outfile << std::setprecision(12);

	// extract the folder of this nvm file
	string nvmfolder(stlplus::folder_part(path));

	// header
	outfile << "NVM_V3" << endl;

	for (auto model : models)
		outfile << model;

	// empty model at the end
	outfile << "0";
	outfile.close();

	LOG(INFO)<< "saved " << models.size() << " models to <"<< path << ">";

}

void NVMReader::saveAsPly(const NVM_Model& model, const char* file) {
	ofstream ofstr;
	ofstr.open(file);
	ofstr << "ply" << '\n' << "format ascii 1.0" << '\n' << "element vertex "
			<< (int) model.points.size() << '\n' << "property float x" << '\n' << "property float y"
			<< '\n' << "property float z" << '\n' << "property uchar diffuse_red" << '\n'
			<< "property uchar diffuse_green" << '\n' << "property uchar diffuse_blue" << '\n'
			<< "end_header" << '\n';

	for (const auto& p : model.points) {
		ofstr << p.xyz[0] << " ";
		ofstr << p.xyz[1] << " ";
		ofstr << p.xyz[2] << " ";

		ofstr << p.rgb[0] << " ";
		ofstr << p.rgb[0] << " ";
		ofstr << p.rgb[0] << "\n";
	}

	ofstr.close();
}

} /* namespace mo3d */
