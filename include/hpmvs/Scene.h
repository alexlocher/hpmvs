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

#ifndef HPMVS_SCENE_H_
#define HPMVS_SCENE_H_

#include <vector>
#include <map>
#include <mutex>
#include <chrono>

#include <Eigen/Dense>

#include <hpmvs/Image.h>
#include <hpmvs/Camera.h>
#include <hpmvs/Patch3d.h>
#include <hpmvs/doctree.h>
#include <hpmvs/HpmvsOptions.h>
#include <hpmvs/NVMReader.h>


namespace mo3d {

class Scene {
public:
	Scene();
	virtual ~Scene();

	bool addCameras(const NVM_Model& model, const HpmvsOptions& options);

	bool initPatches(const NVM_Model& model, const HpmvsOptions& options);

	bool extractCoVisiblilty(const NVM_Model& model, const HpmvsOptions& options);


	Eigen::Vector3f getColor(const Patch3d& patch) const;
	inline Eigen::Vector3f getColor(const int index, const float fx, const float fy, const int level) const;
//	inline Vec3f getColor_old(const int index, const float fx, const float fy, const int level) const;



	int getLevelSupport(const Patch3d& patch, const int minLevel);



	//bool loadFromNVM(const std::string& nvm_fn, const HpmvsOptions& options);

	// fields
	std::map<std::string, int> dict_; // <img name <=> internal id>
	std::vector<Camera> cameras_;
	std::vector<Image> images_;
	std::vector<std::vector<int> > covis_; //

	// --------------------------------------------------------------------------------
	// depth images [img index, level, data ]
	std::vector<std::vector<std::unique_ptr<Eigen::MatrixXf> > > m_depths;
	std::vector<std::vector<std::unique_ptr<std::mutex> > > m_depth_lock;
	static float MAX_DEPTH;
	const double DEPTH_SUBSAMPLE = 2;

	void setDepths(const Patch3d& patch, bool subtract = false);
	float getDetphAtLevel(const int imgIdx, const int xx, const int yy, const int level = 0,
	bool lock = false);
	float getFullDepth(const int imgIdx, const int xx, const int yy, bool lock = false);
	void visualizeDepths(const char* folder);

	void saveAsNVM(const char* folder);
	void savePMats(const char* file);
	void savePoseMats(const char* file);
	// --------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------
	// Test functions
	// --------------------------------------------------------------------------------
	int depthTests(const Patch3d& patch, const float margin);

	int viewBlockTest(const Patch3d& patch, const float margin);

	bool depthTest(const Patch3d& patch, const int image, const float margin, bool neighbours =
	false, bool viewBlock = false);

	bool depthTest(const Patch3d& patch, const int ix, const int iy, const float depth,
			const int image, const float margin, bool viewBlock);

	int pixelFreeTests(const Patch3d& patch);

	bool pixelFreeTest(const Patch3d& patch, const int image);

	// --------------------------------------------------------------------------------

	DynOctTree<Ppatch3d> patchTree_;


	  //----------------------------------------------------------------------
	  // Timing stuff
	  //----------------------------------------------------------------------
	  std::chrono::time_point<std::chrono::steady_clock> t_start; // start point for timing measurements
	  double t_offtime; // duration spent on tasks not relevant for timing
};



Eigen::Vector3f Scene::getColor(const int index, const float fx, const float fy, const int level) const {
	return images_[index].getColor(fx,fy,level);
};


} /* namespace mo3d */

#endif /* HPMVS_SCENE_H_ */
