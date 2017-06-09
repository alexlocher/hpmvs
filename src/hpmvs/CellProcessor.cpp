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

#include <hpmvs/CellProcessor.h>
#include <glog/logging.h>
#include <Eigen/Dense>
#include <set>
#include <stlplus3/file_system.hpp>
#include <iomanip>

#include <hpmvs/Scene.h>

using namespace std;
using namespace Eigen;

namespace mo3d {

CellProcessor::CellProcessor(Scene* scene, const HpmvsOptions& options) :
		f_patchevent(), scene_p(scene), options_p(&options), optimizer_p(0), ptree(0), borderCellFn_(
				0) {
}

CellProcessor::~CellProcessor() {
}

void CellProcessor::filter(Leaf<Ppatch3d>* cell) {
	const int nrPatches = cell->data.size();
	if (nrPatches <= 1)
		return;

	// if we have multiple patches in this node,
	// keep the one with best support among each other
	pair<float, Ppatch3d> best(std::numeric_limits<float>::max(), Ppatch3d());

	for (int ii = 0; ii < nrPatches; ii++) {
		float dist = 0;
		mo3d::Ppatch3d patchA = cell->data[ii];
		Eigen::Vector3f n = patchA->normal_.head(3);
		Eigen::Vector3f x0 = patchA->center_.head(3);
		n.normalize();
		for (int jj = 0; jj < nrPatches; jj++) {
			if (ii == jj)
				continue;
			mo3d::Ppatch3d patchB = cell->data[jj];
			dist += n.dot(patchB->center_.head(3) - x0);
		}
		dist /= (nrPatches - 1);
		if (dist < best.first) {
			best = std::make_pair(dist, patchA);
		}
	}

	// make sure we handle the depths
	for (const Ppatch3d& p : cell->data) {
		if (p != best.second) {
			scene_p->setDepths(*p, true);
			p->images_.clear();
			if (f_patchevent)
				f_patchevent(p);
		}
	}

	cell->data.clear();
	cell->data.push_back(best.second);
}

void CellProcessor::extend(Leaf<Ppatch3d>* cell) {
	if (cell->empty())
		return;
	// check the neighboring cells
	mo3d::Ppatch3d& p = cell->data[0];
	if (p->expanded_) {
		LOG(WARNING)<< "cell is already expanded";
		return;
	}

	// reference image index
	int refImg = p->images_[0];

	// get x and y axis
	Eigen::Vector3f n = p->normal_.head(3);
	Eigen::Vector3f imgX = scene_p->cameras_[refImg].xAxis_;
	Vector3f yaxis = n.cross(imgX).normalized();
	Vector3f xaxis = yaxis.cross(n);

	// extend the patch into N directions
	const int N = 6;
	float extend = cell->width_;
	for (int ii = 0; ii < N; ii++) {
		const float angle = 2.0 * M_PI / N * ii;
		const float dx = cos(angle);
		const float dy = sin(angle);

		// create a new patch
		Ppatch3d newP(new mo3d::Patch3d);
		*newP = *p;
		newP->center_[0] = p->center_[0] + (dx * xaxis[0] + dy * yaxis[0]) * extend;
		newP->center_[1] = p->center_[1] + (dx * xaxis[1] + dy * yaxis[1]) * extend;
		newP->center_[2] = p->center_[2] + (dx * xaxis[2] + dy * yaxis[2]) * extend;
		newP->scale_3dx_ = cell->width_ * 0.9 / 2.0;
		newP->expanded_ = false;
		newP->flatness_ = -1.0;

		// check the corresponding cell
		bool inSide = ptree->getRoot()->contains(newP->center_.head(3));
		Leaf<Ppatch3d>* newLeaf = ptree->at(newP->x(), newP->y(), newP->z());
		if (inSide && (!newLeaf->empty() || newLeaf->width_ < cell->width_)) {
			continue;
		}

		// optimize it
		bool patchGood = optimizer_p->optimize(*newP);
		patchGood = patchGood && newP->scale_3dx_ * 2.0 < cell->width_
				&& newP->scale_3dx_ * 2.0 > cell->width_ / 2.0;
		patchGood = patchGood && (newP->center_ - p->center_).norm() < cell->width_ * 1.5; // should not drift to far

		patchGood = patchGood
				&& scene_p->depthTests(*newP, options_p->DEPTH_TEST_FACTOR)
						>= options_p->MIN_IMAGES_PER_PATCH;
		patchGood = patchGood
				&& scene_p->viewBlockTest(*newP, options_p->DEPTH_TEST_FACTOR)
						< options_p->MIN_IMAGES_PER_PATCH;
		const int nrImgsWithFreePixel = scene_p->pixelFreeTests(*newP);
		patchGood = patchGood && nrImgsWithFreePixel >= options_p->MIN_IMAGES_PER_PATCH - 1
				&& nrImgsWithFreePixel * 1.0 / newP->images_.size() > 0.75;

		if (!patchGood)
			continue;

		if (ptree->getRoot()->contains(newP->center_.head(3)) == false) {
			// this is actually a border cell => handle it
			if (borderCellFn_) {
				float newPriority = (ptree->nodeLevel(cell) + cell->data[0]->priorityReduction_)
						* 10.0;
				(*borderCellFn_)(newP, newPriority);
			}
		} else if (ptree->addConditional(newP, cell->width_ * 0.9, &newLeaf)) {

			if (ptree->nodeLevel(newLeaf) != ptree->nodeLevel(cell)) {
				VLOG(2) << "cell level: " << ptree->nodeLevel(cell) << " ; extended level %d"
									<< ptree->nodeLevel(newLeaf);
				VLOG(2) << "newLeaf Width = " << cell->width_ * 0.75 << " cell width = "
									<< newLeaf->width_;
				exit(1);
			}

			// depth
			scene_p->setDepths(*newP, false);

			if (f_patchevent)
				f_patchevent(newP);

			// add the cell to the processing queue
			processing_queue.push(
					make_pair(
							(ptree->nodeLevel(newLeaf) + cell->data[0]->priorityReduction_) * 10.0,
							newLeaf));
		}
	}
	p->expanded_ = true;
}

bool CellProcessor::refine(Leaf<Ppatch3d>* cell) {
	if (cell->empty())
		return false;
	// check the neighboring cells
	Ppatch3d& p = cell->data[0];

	// reference image index
	int refImg = p->images_[0];

	// reduce patch size
	p->scale_3dx_ = cell->width_ * 0.45 / 2.0;

	// optimize it
	bool patchGood = optimizer_p->optimize(*p);
	patchGood = patchGood && p->scale_3dx_ * 2.0 < cell->width_ / 2.0
			&& p->scale_3dx_ * 2.0 > cell->width_ / 2.0 / 2.0;
	patchGood = patchGood && (ptree->at(p->x(), p->y(), p->z()) == cell);
	patchGood = patchGood
			&& scene_p->depthTests(*p, options_p->DEPTH_TEST_FACTOR)
					>= options_p->MIN_IMAGES_PER_PATCH
			&& scene_p->viewBlockTest(*p, options_p->DEPTH_TEST_FACTOR)
					< options_p->MIN_IMAGES_PER_PATCH;

	if (!patchGood) {
		scene_p->setDepths(*cell->data[0], true);
		ptree->remove(cell);
	}
	return patchGood;
}

void CellProcessor::branch(Leaf<Ppatch3d>* cell) {
	if (cell->empty())
		return;
	// check the neighboring cells
	mo3d::Ppatch3d& p = cell->data[0];
	CHECK (p->expanded_) << "Cell is not expanded => should not branch";

	// reference image index
	int refImg = p->images_[0];

	// get the level support of the patch
	if (scene_p->getLevelSupport(*p, options_p->MINLEVEL) < 1) {
		// consider this cell as exhausted
		return;
	}

	// get x and y axis
	Eigen::Vector3f n = p->normal_.head(3);
	Eigen::Vector3f imgX = scene_p->cameras_[refImg].xAxis_;
	Vector3f yaxis = n.cross(imgX).normalized();
	Vector3f xaxis = yaxis.cross(n);

	// branch the patch into N directions
	const int N = 4;
	float extend = cell->width_ / 4.0;
	vector<Ppatch3d> newPatches;
	for (int ii = 0; ii < N; ii++) {
		const float angle = 2.0 * M_PI / N * ii + M_PI / 4;
		const float dx = cos(angle);
		const float dy = sin(angle);

		// create a new patch
		Ppatch3d newP(new mo3d::Patch3d);
		*newP = *p;
		newP->center_[0] = p->center_[0] + (dx * xaxis[0] + dy * yaxis[0]) * extend;
		newP->center_[1] = p->center_[1] + (dx * xaxis[1] + dy * yaxis[1]) * extend;
		newP->center_[2] = p->center_[2] + (dx * xaxis[2] + dy * yaxis[2]) * extend;
		newP->scale_3dx_ = cell->width_ * 0.45 / 2.0;
		newP->expanded_ = false;
		newP->flatness_ = -1.0;

		// only branch to points staying within the current cell
		if (!cell->contains(Vector3f(newP->x(), newP->y(), newP->z())))
			continue;

		// optimize it
		bool patchGood = optimizer_p->optimize(*newP);

		if (!patchGood)
			continue;

		// check again that the patch is still within the cell (might have moved)
		if (!cell->contains(Vector3f(newP->x(), newP->y(), newP->z())))
			continue;

		// add the patch
		newPatches.push_back(newP);
	}

	// if we did not successfully branched into at least one new patch, consider
	// this cell as exhausted and keep the existing patch
	if (ptree->nodeLevel(cell) >= options_p->PATCH_FINAL_MINLEVEL  && newPatches.size() == 0)
		return;

	// now change the cell
	vector<Ppatch3d> oldPatch;
	Branch<Ppatch3d>* newBranch = cell->split(oldPatch); // from here on, the cell* is not valid
	for (auto oldP : oldPatch) {
		scene_p->setDepths(*oldP, true);
		if (f_patchevent) {
			oldP->images_.clear();
			f_patchevent(oldP);
		}
	}
	oldPatch.clear();

	//
	delete cell;
	cell = nullptr;

	// insert the new patches
	std::set<Leaf<Ppatch3d>*> newCells;
	for (auto& p : newPatches) {
		Leaf<Ppatch3d>* leaf = newBranch->at(p->x(), p->y(), p->z());
		leaf->data.push_back(p);
		scene_p->setDepths(*p, false);
		newCells.insert(leaf);
		if (f_patchevent)
			f_patchevent(p);
	}

	// add the cell to the processing queue
	for (auto it = newCells.begin(); it != newCells.end(); it++) {
		processing_queue.push(
				make_pair((ptree->nodeLevel(*it) + (*it)->data[0]->priorityReduction_) * 10.0,
						*it));
	}
}

void CellProcessor::regularize(Leaf<Ppatch3d>* cell) {
	if (cell->empty())
		return;
	// check the neighboring cells
	Ppatch3d& p = cell->data[0];
	if (!p->expanded_)
		return;

	// reference image index
	int refImg = p->images_[0];

	// get x and y axis
	Eigen::Vector3f n = p->normal_.head(3);
	Eigen::Vector3f imgX = scene_p->cameras_[refImg].xAxis_;
	Vector3f yaxis = n.cross(imgX).normalized();
	Vector3f xaxis = yaxis.cross(n);

	// extend the patch into 4 directions
	float hwin = 2;
	Vector3f c(p->x(), p->y(), p->z());
	std::set<Leaf<Ppatch3d>*> neighCells;
	for (int yy = -hwin; yy <= hwin; yy++) {
		for (int xx = -hwin; xx <= hwin; xx++) {
//				for (int zz = -hwin; zz <= hwin; zz++){
			if (xx == 0 && yy == 0 /* && zz == 0 */)
				continue;
			// calculate the 3d pos of a possible patch extension
			Vector3f c_ext = c + (xx * xaxis + yy * yaxis) * cell->width_;
//					Vector3f c_ext = cell->c_ + Vector3f(xx* cell->width_, yy*cell->width_, zz* cell->width_);
			Leaf<Ppatch3d>* possibleLeaf = ptree->at(c_ext[0], c_ext[1], c_ext[2]);
			if (!possibleLeaf->empty())
				neighCells.insert(possibleLeaf);
//				}
		}
	}

	int nrCells = neighCells.size();
	if (nrCells < 1) {
		p->flatness_ = 2.6;
		return;
	} else if (nrCells < 4) {
		p->flatness_ = 2.5;
		return;
	}

	// calculate the RMS error of this patch
	float dist = 0;
	Eigen::Vector3f x0 = p->center_.head(3);
	n.normalize();
	for (Leaf<Ppatch3d>* nCell : neighCells) {
		Ppatch3d patchB = nCell->data[0];
//		Vector3f nOther(patchB->m_normal[0], patchB->m_normal[1], patchB->m_normal[2]);
		float error = n.dot(patchB->center_.head(3) - x0);
		dist += error * error;
	}
	dist = std::sqrt(dist / nrCells) / cell->width_;

	p->flatness_ = dist;
}

void CellProcessor::processCell(Leaf<Ppatch3d>* cell, float priority) {
	if (cell->data.empty())
		return; // no need for processing

	// stop processing
	if (priority >= (options_p->MAX_TREE_LEVEL + 1) * 10)
		return;

	if (cell->data.size() > 1)
		filter(cell);

	if (!cell->data[0]->expanded_) {
		extend(cell);

		// re-add the cell to the queue (with lower priority) such that it regularized later
		processing_queue.push(
				make_pair((ptree->nodeLevel(cell) + cell->data[0]->priorityReduction_) * 10.0 + 1,
						cell));
		return;
	}

	float flatness = cell->data[0]->flatness_;
	if (flatness < 0) {
		regularize(cell);
//		cell->data[0]->priorityReduction_ += (cell->data[0]->flatness_ < (ptree->nodeLevel(cell)/60.0)) ? 1 : 0;
//		cell->data[0]->priorityReduction_ = (cell->c_ - Vector3f(-3.16, 0.47, 4.79)).norm() > 1.0 ? 2 : 0;
//		cell->data[0]->priorityReduction_ = (cell->c_ - Vector3f(8.28, -4.55, 4.63)).norm() > 3.0 ? 2 : 0;
//		if (ptree->nodeLevel(cell) > 9)
//			cell->data[0]->priorityReduction_ += (cell->c_ - Vector3f(3.45, 0.68, 2.75)).norm() > 0.8 ? 3 : 0; // head der hass
//		else
		cell->data[0]->priorityReduction_ = 0;
		float newPrio = (ptree->nodeLevel(cell) + cell->data[0]->priorityReduction_) * 10.0 + 2;
//		if (newPrio >= (maxTreeLevel*10))
//			return;
		processing_queue.push(make_pair(newPrio, cell));
		return;
	}

	if (flatness > 2.4) {
		scene_p->setDepths(*cell->data[0], true);
		cell->data[0]->images_.clear();
		if (f_patchevent)
			f_patchevent(cell->data[0]);
		ptree->remove(cell);
	} else {
//		if (priority < (maxTreeLevel*10))
		branch(cell);
//		else
//			refine(cell);
	}

}

void CellProcessor::initFromTree(DynOctTree<Ppatch3d>* tree,
		std::function<void(Ppatch3d, const float)>* borderPatch, bool skip_clean) {
	// not used anymore:
	const int minLevel = -1;
	const int maxLevel = std::numeric_limits<int>::max();

	this->ptree = tree;
	this->borderCellFn_ = borderPatch;

	// just to be sure
	while (!processing_queue.empty())
		processing_queue.pop();

	// fill the cells into the queue
	int total = 0, skipped = 0;
	Leaf_iterator<Ppatch3d> end = tree->end();
	for (Leaf_iterator<Ppatch3d> leaf = tree->begin(); leaf != end; leaf++) {
		int nodeLevel = tree->nodeLevel(&(*leaf));
		if (nodeLevel <= maxLevel && nodeLevel >= minLevel && !leaf->empty()) {
			total++;
			// check if this cell is already preprocessed
			if (skip_clean && leaf->data.size() == 1 && leaf->data[0]->dirty_ == false && leaf->data[0]->expanded_ == true){
				skipped++;continue;
			}

			processing_queue.push(std::make_pair(nodeLevel * 10, &(*leaf)));
			if (f_patchevent) {
				for (auto& p : leaf->data)
					f_patchevent(p);
			}
		}
	}

}

bool CellProcessor::processQueue(PatchOptimizer* optimizer, float maxPriority) {
	optimizer_p = optimizer;

	if (ptree == 0)
		return false;

	// handle border patches from the round before...
	bool borderCellsAdded = processBorderCellQueue();

	if (processing_queue.empty())
		return borderCellsAdded;

	float currentPriority = processing_queue.top().first;
	int processedCells = 0;

	// process queue until stopping criteria is met
	while (!processing_queue.empty() && currentPriority <= maxPriority) {
		currentPriority = processing_queue.top().first;
		Leaf<mo3d::Ppatch3d>* newCell = processing_queue.top().second;
		processing_queue.pop();

		// do the processing
		processCell(newCell, currentPriority);

		processedCells++;
	}

	return borderCellsAdded || processedCells > 0;
}

bool CellProcessor::insertBorderCell(Ppatch3d& patch, const float prioritiy) {
	// check if the patch fits the current tree
	if (ptree == nullptr)
		return false;
	if (!ptree->getRoot()->contains(patch->center_.head(3)))
		return false;

	// alright, insert into the queue
	std::lock_guard<std::mutex> queue_lock(borderCellQueueMtx_);
	borderCellQueue_.push(std::pair<float, Ppatch3d>(prioritiy, patch));
	return true;
}

bool CellProcessor::processBorderCellQueue() {
	std::unique_lock<std::mutex> borderQueueLock(borderCellQueueMtx_);
	int nrBorderCellsAdded = 0;
	int nrBorderCells = borderCellQueue_.size();
	while (!borderCellQueue_.empty()) {
		float prio = borderCellQueue_.front().first;
		Ppatch3d newP = borderCellQueue_.front().second;
		borderCellQueue_.pop();

		// the patch is already optimized and everything => now we just have to add it to the tree
		Leaf<Ppatch3d>* newLeaf = ptree->at(newP->x(), newP->y(), newP->z());
		if (ptree->addConditional(newP, newP->scale_3dx_ * 2.0, &newLeaf)) {

			// prevent regularization on border cells
			newP->flatness_ = 0;

			// depth
			scene_p->setDepths(*newP, false);

			if (f_patchevent)
				f_patchevent(newP);

			// add the cell to the processing queue
			processing_queue.push(make_pair(prio, newLeaf));

			nrBorderCellsAdded++;
		}

	}
	borderQueueLock.unlock();
	return nrBorderCellsAdded > 0;
}

void CellProcessor::distributeBorderCell(std::vector<std::unique_ptr<CellProcessor>>* processors,
		Ppatch3d patch, const float priority) {
	// just loop through all available processors...
	for (std::unique_ptr<CellProcessor>& p : *processors) {
		if (p->insertBorderCell(patch, priority))
			return;
	}
}

bool CellProcessor::haveWork() const {
	return !borderCellQueue_.empty() || !processing_queue.empty();
}

}
