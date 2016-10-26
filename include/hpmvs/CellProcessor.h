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

#ifndef HPMVS_CELLPROCESSOR_H_
#define HPMVS_CELLPROCESSOR_H_

#include <queue>
#include <mutex>
#include <utility>
#include <functional>

#include <hpmvs/Patch3d.h>
#include <hpmvs/doctree.h>
#include <hpmvs/PatchOptimizer.h>

namespace mo3d {

typedef std::pair<float, Leaf<Ppatch3d>*> PriorityCell;

class CellCompare {
public:
	bool operator()(const PriorityCell& lhs, const PriorityCell& rhs) const {
		if (lhs.first == rhs.first) {
			if (lhs.second->data.size() > 0 && rhs.second->data.size() > 0)
				return lhs.second->data[0]->ncc_ < rhs.second->data[0]->ncc_;
		}
		return lhs.first > rhs.first;
	}
};

class CellProcessor {
private:

	Scene* scene_p;

	// the processing queue
	std::priority_queue<PriorityCell, std::deque<PriorityCell>, CellCompare> processing_queue;

	// callback function for patchchanges
	std::function<void(const mo3d::Ppatch3d&)> f_patchevent;

	// the options
	const HpmvsOptions* options_p;

	// stuff related to the current tree
	DynOctTree<Ppatch3d>* ptree;
	std::function<void(Ppatch3d, const float)>* borderCellFn_;
	PatchOptimizer* optimizer_p;

	std::queue<std::pair<float, Ppatch3d> > borderCellQueue_;
	std::mutex borderCellQueueMtx_;

public:
	CellProcessor(Scene* scene, const HpmvsOptions& options);
	virtual ~CellProcessor();

	void initFromTree(DynOctTree<Ppatch3d>* tree,
			std::function<void(Ppatch3d, const float)>* borderPatch = 0,
 			bool skip_clean = false);

	bool processQueue(PatchOptimizer* optimizer, float maxPriority =
			std::numeric_limits<float>::max());

	bool insertBorderCell(Ppatch3d& patch, const float priority);

	bool haveWork() const;

	static void distributeBorderCell(std::vector<std::unique_ptr<CellProcessor>>* processors,
			Ppatch3d patch, const float priority);

private:

	void filter(Leaf<Ppatch3d>* cell);
	void extend(Leaf<Ppatch3d>* cell);
	bool refine(Leaf<Ppatch3d>* cell);
	void regularize(Leaf<Ppatch3d>* cell);
	void branch(Leaf<Ppatch3d>* cell);

	void processCell(Leaf<Ppatch3d>* cell, float priority);

	bool processBorderCellQueue();

};

} // namespace

#endif /* HPMVS_CELLPROCESSOR_H_ */
