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

#include <stlplus3/portability/portability.hpp>

#include <hpmvs/NVMReader.h>
#include <hpmvs/CellProcessor.h>
#include <hpmvs/Scene.h>

#include <omp.h>

using namespace std;

// FLAGS
DEFINE_string(nvm, "", "input nvm file");
DEFINE_string(outdir, "/tmp/hpmvs", "output directory");
DEFINE_bool(forcelogtostderr, true, "log to stderr");
DEFINE_int32(subtrees, 100, "min number subtrees the model is split into"	);

template<class Element>
void getSubTrees(DynOctTree<Element>& tree,
		std::vector<std::shared_ptr<DynOctTree<Element> > >& subTrees, const int minTrees = 2) {

	if (minTrees < 2){
		subTrees.emplace_back(std::shared_ptr<DynOctTree<Element> >(&tree));
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
			return;

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
	std::vector<mo3d::NVM_Model> models;
	mo3d::NVMReader::readFile(dataset.c_str(), models, true);
	if (models.empty()) {
		LOG(WARNING)<< "no models found in NVM file";
		return EXIT_FAILURE;
	}

	scene.addCameras(models[0], options);
	scene.extractCoVisiblilty(models[0], options);
	scene.initPatches(models[0], options);
	scene.savePMats(stlplus::create_filespec(options.OUTFOLDER, "pmats", "txt").c_str());
	scene.savePoseMats(stlplus::create_filespec(options.OUTFOLDER, "poses", "txt").c_str());

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
	const int maxPrio = 200;
	for (int prio = 0; prio < maxPrio; prio += 1) {

		std::atomic<uint32_t> trees_changed(0);

#pragma omp parallel for schedule(dynamic)
		for (int ii = 0; ii < subTrees.size(); ii++) {
			const int threadId = omp_get_thread_num();
			if (cellProcessors[ii]->processQueue(&optimizers[threadId], prio))
				trees_changed++;
		}

		if (trees_changed.load() > 0 && ((int) prio) % 10 < 3) {
//			scene.patchTree_.toPly(
//					stlplus::create_filespec(options.OUTFOLDER,
//							"tree-" + std::to_string(prio) + "-").c_str(),
//					false);
			scene.patchTree_.toExtPly(
					stlplus::create_filespec(options.OUTFOLDER, "patches-" + std::to_string(prio),
							"ply").c_str());

			// save depth debug
//			if (((int) prio) / 10 > 7) {
//				string depthFolder = stlplus::create_filespec(outfolder,
//						"afterL" + std::to_string(prio));
//				if (!stlplus::folder_exists(depthFolder))
//					stlplus::folder_create(depthFolder);
//				scene.visualizeDepths(depthFolder.c_str());
//
//				// -----------------------------------------
//				scene.saveAsNVM(depthFolder.c_str());
//			}

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
			stlplus::create_filespec(options.OUTFOLDER, "patches-" + std::to_string(procTime)+ " sec",
					"ply").c_str());

	// -----------------------------------------------------------------
	return EXIT_SUCCESS;
}

int main(int argc, char* argv[]) {
	google::InitGoogleLogging(argv[0]);
	FLAGS_colorlogtostderr = true;
	FLAGS_logtostderr = FLAGS_forcelogtostderr;
	GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

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
	LOG(INFO)<< "dataset      : <" << FLAGS_nvm << ">";
	LOG(INFO)<< "dataset      : <" << FLAGS_outdir << ">";
	LOG(INFO)<< "dataset      : <" << omp_get_max_threads();

	// set the options
	mo3d::HpmvsOptions options;
	options.OUTFOLDER = FLAGS_outdir;

	// launch the actual thing
	return hp_pmvs(FLAGS_nvm, options);
}

