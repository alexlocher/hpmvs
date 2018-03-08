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

#include <iostream>
#include <string>
#include <functional>

#include <glog/logging.h>
#include <gflags/gflags.h>

#include <stlplus3/file_system.hpp>

#include <nvmtools/NVMFile.h>

#include <hpmvs/CellProcessor.h>
#include <hpmvs/Scene.h>

#include <omp.h>

using namespace std;

// FLAGS
DEFINE_string(nvm, "", "input nvm file");
DEFINE_string(outdir, "/tmp/hpmvs", "output directory");
DEFINE_bool(forcelogtostderr, true, "log to stderr");
DEFINE_int32(subtrees, 100, "min number subtrees the model is split into");
DEFINE_int32(maxtreelevel, 20, "maximum level of the octree");
DEFINE_int32(patch_level_final_min, 9, "in case branching stops, min level to keep the lowres patch");
DEFINE_int32(patch_level_init_max, 9, "the max tree level on which patches are initialized");
DEFINE_bool(more_output, false, "save more intermediate pointclouds");
DEFINE_int32(light_output, 0, "also save a leightweight pointcloud as output (provided int is the level of the output cloud e.g. 80)");
DEFINE_bool(only_sphere, false, "only reconstruct points within a sphere around the scene center");

DEFINE_double(ncc1, 0.6 , "patch NCC threshold before optimization (alpha 1)");
DEFINE_double(ncc2, 0.8 , "patch NCC threshold after optimization (alpha 2)");

DEFINE_double(minAngle, 10, "minimum angle between any two images of a patch");
DEFINE_double(maxAngle, 60, "maximum angle between any two images of a patch");

template<class Element>
void getSubTrees(DynOctTree<Element>& tree,
		std::vector<std::shared_ptr<DynOctTree<Element> > >& subTrees, const int minTrees = 2) {

	if (minTrees < 2) {
		// create a subtree from the root
		subTrees.emplace_back(
				std::shared_ptr<DynOctTree<Element> >(new DynOctTree<Element>(tree.getRoot())));
		return;
	}

	// do a first split
	tree.getSubTrees(subTrees);

	while (subTrees.size() < minTrees) {

		// get the subtree with the most leafs
		int maxIndex = -1;
		int maxLeafs = -1;
		for (int ii = 0; ii < subTrees.size(); ii++) {
			int nrLeafs = subTrees[ii]->getRoot()->nrLeafs(); // recursive and slow!
			if (nrLeafs > maxLeafs) {
				maxLeafs = nrLeafs;
				maxIndex = ii;
			}
		}

		// doesn't make sense to split, if we have to few points
		if (maxLeafs < 100)
			break;

		std::shared_ptr<DynOctTree<Element> > maxTree = subTrees[maxIndex];
		std::vector<std::shared_ptr<DynOctTree<Element> >> newSubTrees;
		maxTree->getSubTrees(newSubTrees); // works because we are working with pointers
		for (int ii = 0; ii < subTrees.size(); ii++)
			if (ii != maxIndex)
				newSubTrees.push_back(subTrees[ii]);

		subTrees.swap(newSubTrees);

	}

	LOG(INFO)<< "Split to " << subTrees.size() << " subtrees";
	for (int ii = 0; ii < subTrees.size(); ii++) {
		LOG(INFO)<< "Nr " << ii << " => " << subTrees[ii]->getRoot()->nrLeafs() << " Leafs";
	}
}


int hp_pmvs(const std::string& dataset, const mo3d::HpmvsOptions options) {

	// create a scene, holding everything together
	mo3d::Scene scene;
	std::vector<mo3d::Ppatch3d> initPatches_new;

	// Initialize the scene from the nvm model
	std::vector<nvmtools::NVM_Model> models;
	nvmtools::NVMFile::readFile(dataset.c_str(), models, true);
	if (models.empty()) {
		LOG(WARNING)<< "no models found in NVM file";
		return EXIT_FAILURE;
	}

	scene.addCameras(models[0], options);
	scene.extractCoVisiblilty(models[0], options);
	mo3d::HpmvsOptions initOptions = options;
	initOptions.START_LEVEL = 2;
	scene.initPatches(models[0], options);
	if (FLAGS_more_output)
		scene.patchTree_.toExtPly(stlplus::create_filespec(options.OUTFOLDER, "patches-init" ,"ply").c_str());

	const int nrThreads = omp_get_max_threads();

	// create optimizers
	std::vector<mo3d::PatchOptimizer> optimizers;
	for (int ii = 0; ii < nrThreads; ii++)
		optimizers.emplace_back(options, &scene);

	// divide the initial tree into multiple subtrees
	std::vector<std::shared_ptr<DynOctTree<mo3d::Ppatch3d> > > subTrees;
	getSubTrees(scene.patchTree_, subTrees, FLAGS_subtrees);

	// initialize the CellProcessors
	std::vector<std::unique_ptr<mo3d::CellProcessor> > cellProcessors;
	std::function<void(mo3d::Ppatch3d, const float)> borderPatchFn = std::bind(
			&mo3d::CellProcessor::distributeBorderCell, &cellProcessors, std::placeholders::_1,
			std::placeholders::_2);
	for (int ii = 0; ii < subTrees.size(); ii++) {
		cellProcessors.emplace_back(new mo3d::CellProcessor(&scene, options));
		cellProcessors.back()->initFromTree(subTrees[ii].get(), &borderPatchFn);
	}

	// some timing
	const auto start = std::chrono::steady_clock::now();

	// process the queues level by level
	const int maxPrio = (options.MAX_TREE_LEVEL + 1) * 10;
	for (int prio = 0; prio < maxPrio; prio += 1) {

		std::atomic<uint32_t> trees_changed(0);

#pragma omp parallel for schedule(dynamic)
		for (int ii = 0; ii < subTrees.size(); ii++) {
			const int threadId = omp_get_thread_num();
			if (cellProcessors[ii]->processQueue(&optimizers[threadId], prio))
				trees_changed++;
		}

		if (trees_changed.load() > 0 && ((int) prio) % 10 < 3) {
			if (((int) prio) % 10 == 0 || FLAGS_more_output)
				scene.patchTree_.toExtPly(
					stlplus::create_filespec(options.OUTFOLDER, "patches-" + std::to_string(prio),
							"ply").c_str());

			if (((int)prio) == FLAGS_light_output && FLAGS_light_output > 0){
				scene.patchTree_.toExtPly(stlplus::create_filespec(options.OUTFOLDER, "patches-light","ply").c_str(),
						true, // binary
						false, // no normals
						false, // no scale
						false); // no visibility
			}

			LOG(INFO)<< "prio " << prio << " finished";
		}

		// are we finished ?
		bool moreWork = false;
		for (int ii = 0; ii < cellProcessors.size(); ii++)
			moreWork |= cellProcessors[ii]->haveWork();

		if (moreWork == false)
			break;
	}
	scene.patchTree_.cellHistogram();
	const auto end = std::chrono::steady_clock::now();
	double procTime = std::chrono::duration<double>(end - start).count();
	LOG(INFO)<< "Done within " << procTime << " seconds";

	// save the thing:
//	scene.patchTree_.toPly(
//			stlplus::create_filespec(options.OUTFOLDER,
//					"tree-" + std::to_string(procTime) + " sec -").c_str(),
//			false);
	scene.patchTree_.toExtPly(
			stlplus::create_filespec(options.OUTFOLDER,
					"patches-final", "ply").c_str());

	if (FLAGS_light_output > 0){
		scene.patchTree_.toExtPly(stlplus::create_filespec(options.OUTFOLDER, "patches-final-light","ply").c_str(),
				true, // binary
				false, // no normals
				false, // no scale
				false); // no visibility
	}

	// -----------------------------------------------------------------
	return EXIT_SUCCESS;
}

int main(int argc, char* argv[]) {
	google::InitGoogleLogging(argv[0]);
	FLAGS_colorlogtostderr = true;
	FLAGS_logtostderr = FLAGS_forcelogtostderr;
	gflags::ParseCommandLineFlags(&argc, &argv, true);

	LOG(INFO)<< " =============================================================== ";
	LOG(INFO)<< " ======== welcome to the Progressive Multiview Stereo    ======= ";
	LOG(INFO)<< " =============================================================== ";

	// check the arguments
	CHECK(stlplus::file_readable(FLAGS_nvm)) << "input file <" << FLAGS_nvm << "> not readable";
	CHECK(stlplus::folder_exists(FLAGS_outdir) || stlplus::folder_create(FLAGS_outdir))
																									<< "Unable to create output folder <"
																									<< FLAGS_outdir
																									<< ">";

	// output some information
	LOG(INFO)<< "dataset        : <" << FLAGS_nvm << ">";
	LOG(INFO)<< "out directory  : <" << FLAGS_outdir << ">";
	LOG(INFO)<< "number threads : <" << omp_get_max_threads() << ">";

	// set the options
	mo3d::HpmvsOptions options;
	options.OUTFOLDER = FLAGS_outdir;
	options.MAX_TREE_LEVEL = FLAGS_maxtreelevel;
	options.PATCH_FINAL_MINLEVEL = FLAGS_patch_level_final_min;
	options.PATCH_INIT_MAXLEVEL = FLAGS_patch_level_init_max;
	options.FILTER_SCENE_CENTER = FLAGS_only_sphere;
	options.NCC_ALPHA_1 = FLAGS_ncc1;
	options.NCC_ALPHA_2 = FLAGS_ncc2;
	options.MIN_ANGLE = FLAGS_minAngle;
	options.MAX_ANGLE = FLAGS_maxAngle;

	// launch the actual thing
	return hp_pmvs(FLAGS_nvm, options);
}
