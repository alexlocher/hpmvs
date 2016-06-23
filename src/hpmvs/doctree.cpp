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

#include <hpmvs/doctree.h>


// ---------------------------------------------------------------------
// Cell implementation
// ---------------------------------------------------------------------

Cell::Cell(const Eigen::Vector3f c, const float width) : c_(c), width_(width), parent_(nullptr), parentIdx_(0),  type_(UNKNOWN) {};

Cell::Cell(Cell* parent, unsigned int idx) :
		parent_(parent), parentIdx_(idx), type_(UNKNOWN) {
	this->width_ = parent_->width_ / 2.0;
	this->c_[0] = parent_->c_[0] + (!!(idx & (1 << 0)) ? 1.0 : -1.0) * this->width_/2.0;
	this->c_[1] = parent_->c_[1] + (!!(idx & (1 << 1)) ? 1.0 : -1.0) * this->width_/2.0;
	this->c_[2] = parent_->c_[2] + (!!(idx & (1 << 2)) ? 1.0 : -1.0) * this->width_/2.0;
}

bool Cell::contains(const Eigen::Vector3f& p) const {
	const float hw = width_/2.0;
	return  p[0] > c_[0] - hw && p[1] > c_[1] - hw && p[2] > c_[2] - hw &&
			p[0] <= c_[0] + hw && p[1] <= c_[1] + hw && p[2] <= c_[2] + hw;
}
