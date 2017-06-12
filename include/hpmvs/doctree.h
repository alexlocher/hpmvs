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

#ifndef HPMVS_DOCTREE_H
#define HPMVS_DOCTREE_H

#include <Eigen/Dense>
#include <glog/logging.h>
#include <fstream>
#include <list>
#include <memory>

template<typename Element> class Leaf;
template<typename Element> class Branch;
template<typename Element> class DynOctTree;

template<typename Element> void getBoundingBox(const std::vector<Element>& e, Eigen::Vector3f& min,
		Eigen::Vector3f& max);

class Cell {
private:
	Cell(const Cell& other) = delete; // non construction-copyable
	Cell& operator=(const Cell&) = delete; // non copyable
protected:
	Cell(const Eigen::Vector3f c, const float width);
	Cell(Cell* parent, unsigned int idx);
public:
	virtual ~Cell() {
	}

	bool contains(const Eigen::Vector3f& p) const;

	enum TYPE {
		UNKNOWN, BRANCH, LEAF,
	};
	TYPE type_;
	Cell* parent_;
	unsigned char parentIdx_;

	Eigen::Vector3f c_; // cells center point
	float width_;		// cells width
};

template<typename Element>
class Branch: public Cell {
private:
	Branch(const Branch& other) = delete; // non construction-copyable
	Branch& operator=(const Branch&) = delete; // non copyable

public:
	Branch(const Eigen::Vector3f c, float width);
	Branch(Cell* parent, unsigned int idx);
	virtual ~Branch();

	Cell* children[8];
//	void collapse(); // recursively collapses all leaves into the cell
//	Cell* collapseIfEmtpy(); // collapses this branch if all children are empty leafs
	bool empty() const;

	size_t nrLeafs() const;

	Leaf<Element>* at(float x, float y, float z); // get a pointer to the deepest value of given point coordinate
};

template<typename Element>
class Leaf: public Cell {
private:
	Leaf(const Leaf& other) = delete; // non construction-copyable
	Leaf& operator=(const Leaf&) = delete; // non copyable
public:
	Leaf(Cell* parent, unsigned int idx);
	virtual ~Leaf();
	bool empty() const {
		return data.empty();
	}
	std::vector<Element> data;
	Branch<Element>* split(std::vector<Element>& data); // does a split
};

template<typename Element>
class Leaf_iterator: public std::iterator<std::forward_iterator_tag, // type of iterator
		Leaf<Element>, ptrdiff_t, Leaf<Element>*, Leaf<Element>&> // Info about iterator
{
private:
	const Branch<Element>* root_;
	Leaf<Element>* node_;

public:
	Leaf_iterator(const Branch<Element>* tree, Leaf<Element>* node);
	Leaf<Element>& operator*();
	Leaf<Element>* operator->();
	Leaf_iterator operator++();
	Leaf_iterator& operator++(int);
	bool equal(Leaf_iterator const& rhs) const;
	bool operator==(const Leaf_iterator<Element>& rhs) {
		return this->equal(rhs);
	}
	bool operator!=(const Leaf_iterator<Element>& rhs) {
		return !this->equal(rhs);
	}

};

template<typename Element>
class DynOctTree {
public:

	DynOctTree();
	DynOctTree(const Eigen::Vector3f c, const float width);
	DynOctTree(Branch<Element>* root, bool isSubTree = true, const int rootLevel = 0);
	~DynOctTree();

	// get a pointer to the deepest value of given point coordinate
	inline Leaf<Element>* at(const float x, const float y, const float z) {
		return root->at(x, y, z);
	}

	inline Branch<Element>* getRoot() {
		return root;
	}

	Leaf_iterator<Element> begin();
	Leaf_iterator<Element> end();

	Leaf<Element>* add(Element e);
	Leaf<Element>* add(Element e, const float width);
	Leaf<Element>* addSingleCell(Element e);
	bool addConditional(Element e, const float width, Leaf<Element>** leaf);
	Leaf<Element>* remove(Element e);
	Leaf<Element>* remove(Leaf<Element>* leaf);

	Branch<Element>* swapRoot(Branch<Element>* newRoot);

	inline int nodeLevel(const Cell* node) const;

	void toPly(const char*, bool noScalar = true);
	void toExtPly(const char*, bool binary = false,	bool normal = true,	bool scale = true, bool visibility = true);

	std::list<int> pathToRoot(const Cell* node);
	std::vector<int> cellHistogram();

	void getSubTrees(std::vector<std::shared_ptr<DynOctTree<Element> > >& trees);

private:
	Leaf<Element>* remove_internal(Leaf<Element>* leaf);

	Branch<Element>* root;
	bool ownsRoot;
	int rootLevel_; // in case of subtree
};

// ---------------------------------------------------------------------
// Leaf implementation
// ---------------------------------------------------------------------
template<typename Element>
Leaf<Element>::Leaf(Cell* parent, unsigned int idx) :
		Cell(parent, idx) {
	this->type_ = Cell::LEAF;
}

template<typename Element>
Leaf<Element>::~Leaf() {
}

template<typename Element>
Branch<Element>* Leaf<Element>::split(std::vector<Element>& d) {
	Branch<Element> *new_branch = new Branch<Element>(this->parent_, this->parentIdx_);
	reinterpret_cast<Branch<Element>*>(this->parent_)->children[this->parentIdx_] = new_branch; // attach the new branch to parent

	d.swap(data);
	this->parent_ = nullptr;
	return new_branch;
}

// ---------------------------------------------------------------------
// Branch implementation
// ---------------------------------------------------------------------
template<typename Element>
Branch<Element>::Branch(Cell* parent, unsigned int idx) :
		Cell(parent, idx) {
	this->type_ = Cell::BRANCH;

	// create children
	for (int ii = 0; ii < 8; ii++)
		this->children[ii] = new Leaf<Element>(this, ii);
}

template<typename Element>
Branch<Element>::Branch(const Eigen::Vector3f c, float width) :
		Cell(c, width) {
	this->type_ = Cell::BRANCH;

	// create children
	for (int ii = 0; ii < 8; ii++)
		this->children[ii] = new Leaf<Element>(this, ii);
}

template<typename Element>
Branch<Element>::~Branch() {
	// recursively delete children
	for (int ii = 0; ii < 8; ii++)
		if (nullptr != this->children[ii])
			delete children[ii];
}

template<typename Element>
bool Branch<Element>::empty() const {
	for (int ii = 0; ii < 8; ii++) {
		if (children[ii]->type_ == Cell::LEAF
				&& !reinterpret_cast<Leaf<Element>*>(children[ii])->empty())
			return false;
		if (children[ii]->type_ == Cell::BRANCH
				&& !reinterpret_cast<Branch<Element>*>(children[ii])->empty())
			return false;
	}
	return true;
}

template<typename Element>
size_t Branch<Element>::nrLeafs() const {
	size_t leafs = 0;
	for (int ii = 0; ii < 8; ii++) {
		if (children[ii]->type_ == Cell::LEAF
				&& !reinterpret_cast<Leaf<Element>*>(children[ii])->empty())
			leafs++;
		else if (children[ii]->type_ == Cell::BRANCH)
			leafs += reinterpret_cast<Branch<Element>*>(children[ii])->nrLeafs();
	}
	return leafs;
}

template<typename Element>
Leaf<Element>* Branch<Element>::at(const float x, const float y, const float z) {
	unsigned char idx = ((z > this->c_[2]) << 2) | ((y > this->c_[1]) << 1) | (x > this->c_[0]);
	if (children[idx]->type_ == Cell::LEAF)
		return reinterpret_cast<Leaf<Element>*>(children[idx]);
	return reinterpret_cast<Branch<Element>*>(children[idx])->at(x, y, z);
}

// ---------------------------------------------------------------------
// Iterator implementation
// ---------------------------------------------------------------------
template<typename Element>
Leaf_iterator<Element>::Leaf_iterator(const Branch<Element>* tree, Leaf<Element>* node) :
		node_(node), root_(tree) {
}

template<typename Element>
Leaf<Element>& Leaf_iterator<Element>::operator*() {
	return *node_;
}

template<typename Element>
Leaf<Element>* Leaf_iterator<Element>::operator->() {
	return node_;
}

template<typename Element>
Leaf_iterator<Element> Leaf_iterator<Element>::operator++() {
	Leaf_iterator<Element> before = *this;
	(*this)++;
	return before;
}

template<typename Element>
Leaf_iterator<Element>& Leaf_iterator<Element>::operator++(int) {
	const Leaf<Element>* before = node_;
	Cell* b = node_->parent_;
	unsigned int idx = node_->parentIdx_;

	// travel up the tree until we can increment the index by one (or we reach root)
	while (idx + 1 > 7 && b != root_) {
		idx = b->parentIdx_;
		b = b->parent_;
	}

	// if we came from idx 7 and are on root now -> this was the last node
	if (b == root_ && idx == 7){
		node_ = nullptr;
		return *this;
	}

	// get the next index (can be branch or leaf)
	idx = (idx > 5) ? 7 : idx + 1; // also handles the root case
	b = reinterpret_cast<Branch<Element>*>(b)->children[idx];

	// possibly travel down the tree until the first leaf
	while (b->type_ != Cell::LEAF) {
		b = reinterpret_cast<Branch<Element>*>(b)->children[0];
	}

	node_ = reinterpret_cast<Leaf<Element>*>(b);

	return *this;
}

template<typename Element>
bool Leaf_iterator<Element>::equal(Leaf_iterator const& rhs) const {
	return (this->root_ == rhs.root_ && this->node_ == rhs.node_);
}

// ---------------------------------------------------------------------
// Octree implementation
// ---------------------------------------------------------------------
template<typename Element>
DynOctTree<Element>::DynOctTree() {
	root = new Branch<Element>(Eigen::Vector3f::Zero(), 1.0);
	ownsRoot = true;
	rootLevel_ = 0;
}

template<typename Element>
DynOctTree<Element>::DynOctTree(const Eigen::Vector3f c, const float width) {
	root = new Branch<Element>(c, width);
	ownsRoot = true;
	rootLevel_ = 0;
	LOG(INFO)<< "initialized tree with center [" << c.transpose() << "] and width = " << width;
}

template<typename Element>
DynOctTree<Element>::DynOctTree(Branch<Element>* root, bool isSubTree, const int rootLevel) {
	this->root = root;
	ownsRoot = !isSubTree;
	rootLevel_ = rootLevel;
}

template<typename Element>
DynOctTree<Element>::~DynOctTree() {
	if (ownsRoot)
		delete root;
}

template<typename Element>
Leaf<Element>* DynOctTree<Element>::add(Element e) {
	root->at(e->x(), e->y(), e->z())->data.push_back(e);
}

template<typename Element>
Leaf<Element>* DynOctTree<Element>::addSingleCell(Element e) {
	Leaf<Element>* leaf = root->at(e->x(), e->y(), e->z());

	if (leaf->empty())
		leaf->data.push_back(e);
	else {
		std::vector<Element> existing;
		Branch<Element>* newBranch = leaf->split(existing);
		existing.push_back(e);
		delete leaf;

		for (Element& els : existing)
			newBranch->at(els->x(), els->y(), els->z())->data.push_back(els);
	}
}

/**
 * adds the element and pushes it as deep as the cells width is smaller than width
 * @param e
 * @param width
 * @return
 */
template<typename Element>
Leaf<Element>* DynOctTree<Element>::add(Element e, const float width) {
	Leaf<Element>* leaf = root->at(e->x(), e->y(), e->z());
	Branch<Element>* branch = nullptr;
	std::vector<Element> buffer;
	while (leaf->width_ / 2.0 > width) {
		branch = leaf->split(buffer);
		delete leaf;
		for (int ii = 0; ii < buffer.size(); ii++)
			branch->at(buffer[ii]->x(), buffer[ii]->y(), buffer[ii]->z())->data.push_back(
					buffer[ii]);
		buffer.clear();
		leaf = branch->at(e->x(), e->y(), e->z());
	}
	leaf->data.push_back(e);
	return leaf;
}

template<typename Element>
bool DynOctTree<Element>::addConditional(Element e, const float width, Leaf<Element>** outleaf) {

	// check if the leaf is empty and not already on a lower level
	Leaf<Element>* leaf = root->at(e->x(), e->y(), e->z());
	if (!leaf->empty() || leaf->width_ < width) {
		*outleaf = leaf;
		return false;
	}

	// add the element to the leaf
	Branch<Element>* branch = nullptr;
	std::vector<Element> buffer;
	while (leaf->width_ / 2.0 > width) {
		branch = leaf->split(buffer);
		delete leaf;
		// we do not have to add any elements from the buffer, because it is empty due to the conditional

		leaf = branch->at(e->x(), e->y(), e->z());
	}
	leaf->data.push_back(e);
	*outleaf = leaf;
	return true;
}

template<typename Element>
Leaf<Element>* DynOctTree<Element>::remove_internal(Leaf<Element>* leaf) {
	Branch<Element>* oldBranch = reinterpret_cast<Branch<Element>*>(leaf->parent_);
	if (leaf->empty() && oldBranch->empty()) {
		Leaf<Element>* newLeaf = new Leaf<Element>(oldBranch->parent_, oldBranch->parentIdx_);
		reinterpret_cast<Branch<Element>*>(oldBranch->parent_)->children[oldBranch->parentIdx_] =
				newLeaf;
		delete oldBranch;
		return newLeaf;
	} else {
		return leaf;
	}
}

template<typename Element>
Leaf<Element>* DynOctTree<Element>::remove(Element e) {
	Leaf<Element>* leaf = root->at(e->x(), e->y(), e->z());
	auto it = std::find(leaf->data.begin(), leaf->data.end(), e);
	if (it == leaf->data.end())
		return leaf;
	leaf->data.erase(it);

	return remove_internal(leaf);
}

template<typename Element>
Leaf<Element>* DynOctTree<Element>::remove(Leaf<Element>* leaf) {
	leaf->data.clear();
	return remove_internal(leaf);
}

template<typename Element>
Branch<Element>* DynOctTree<Element>::swapRoot(Branch<Element>* newRoot) {
	Branch<Element>* oldRoot = root;
	root = newRoot;
	return oldRoot;
}

template<typename Element>
int DynOctTree<Element>::nodeLevel(const Cell* node) const {
	return (int) log2(root->width_ / (node->width_)) + rootLevel_;
}

template<typename Element>
std::list<int> DynOctTree<Element>::pathToRoot(const Cell* node) {
	std::list<int> path;
	const Cell* cell = node;
	while (cell != root) {
		path.push_front(cell->parentIdx_);
		cell = cell->parent_;
	}
	return path;
}

template<typename Element>
Leaf_iterator<Element> DynOctTree<Element>::begin() {
	Branch<Element>* branch = root;
	Cell* node = nullptr;
	while ((node = branch->children[0])->type_ != Cell::LEAF)
		branch = reinterpret_cast<Branch<Element>*>(node);
	return Leaf_iterator<Element>(root, reinterpret_cast<Leaf<Element>*>(node));
}

template<typename Element>
Leaf_iterator<Element> DynOctTree<Element>::end() {
	Branch<Element>* branch = root;
	Cell* node = nullptr;
//	while ((node = branch->children[7])->type_ != Cell::LEAF)
//		branch = reinterpret_cast<Branch<Element>*>(node);
	return Leaf_iterator<Element>(root, reinterpret_cast<Leaf<Element>*>(node));
}

template<typename Element>
std::vector<int> DynOctTree<Element>::cellHistogram() {
	std::vector<int> hist;
	Leaf_iterator<Element> end = this->end();
	for (Leaf_iterator<Element> leaf = this->begin(); leaf != end; leaf++) {
		if (leaf->empty())
			continue;
		int level = nodeLevel(&(*leaf));
		if (level >= hist.size())
			hist.resize(level + 1);
		hist[level]++;
	}

	LOG(INFO)<< "Tree Cell Histogram: ";
	for (int ii = 0; ii < hist.size(); ii++) {
		LOG(INFO)<< "L" << ii << " => " << hist[ii];
	}
	return hist;
}

template<typename Element>
void DynOctTree<Element>::getSubTrees(std::vector<std::shared_ptr<DynOctTree<Element> > >& trees) {
	for (int ii = 0; ii < 8; ii++) {
		if (root->children[ii] != 0 && root->children[ii]->type_ == Cell::BRANCH) {
			Branch<Element>* subRoot = reinterpret_cast<Branch<Element>*>(root->children[ii]);
//			if (!subRoot->empty())
				trees.emplace_back(
						std::make_shared<DynOctTree<Element> >(subRoot, true, rootLevel_ + 1));
		}
	}
}

template<typename Element>
void DynOctTree<Element>::toExtPly(const char* name,bool binary,bool normal,
																		bool scale,	bool visibility) {

	// grab the patches
	std::vector<Element> patches;
	Leaf_iterator<Element> end = this->end();
	for (Leaf_iterator<Element> leaf = this->begin(); leaf != end; leaf++) {
		std::vector<Element>& ps = leaf->data;
		for (const auto& p : ps)
			patches.emplace_back(p);
	}

	// test endianness
	int32_t n;
	bool bigEndian = *(char *) &n == 1;

	// output the header
	std::ofstream pFile(name, std::ofstream::out);
	pFile << "ply" << std::endl;
	if (binary && bigEndian)
		pFile << "format binary_big_endian 1.0" << std::endl;
	else if (binary)
		pFile << "format binary_little_endian 1.0" << std::endl;
	else
		pFile << "format ascii 1.0" << std::endl;
	pFile << "element vertex " << (int) patches.size() << std::endl;
	pFile << "property float x" << std::endl;
	pFile << "property float y" << std::endl;
	pFile << "property float z" << std::endl;
	if (normal){
		pFile << "property float nx" << std::endl;
		pFile << "property float ny" << std::endl;
		pFile << "property float nz" << std::endl;
	}
	pFile << "property uchar red" << std::endl;
	pFile << "property uchar green" << std::endl;
	pFile << "property uchar blue" << std::endl;
	if (scale)
		pFile << "property float scalar_scale" << std::endl;
	if (visibility){
		pFile << "element point_visibility " << (int) patches.size() << std::endl;
		pFile << "property list uint uint visible_cameras" << std::endl;
	}
	pFile << "end_header" << std::endl;
	pFile.close();

	// output the points
	std::ofstream pData(name,
			(binary ? std::ofstream::binary | std::ofstream::app : std::ofstream::app));

	for (const auto& p : patches) {
		if (binary) {
			Eigen::Vector3f v(p->x(), p->y(), p->z());
			pData.write((char*) v.data(), 3 * sizeof(float));

			if (normal){
				v = p->normal_.head(3);
				pData.write((char*) v.data(), 3 * sizeof(float));
			}

			Eigen::Matrix<unsigned char, 3, 1> c(p->color_[0], p->color_[1], p->color_[2]);
			pData.write((char*) c.data(), 3 * sizeof(unsigned char));

			if (scale)
				pData.write((char*) &p->scale_3dx_, sizeof(p->scale_3dx_));

		} else {
			pData << p->x() << " " << p->y() << " " << p->z() << " ";
			if (normal)
				pData << p->normal_[0] << " " << p->normal_[1] << " " << p->normal_[2] << " ";
			pData << (int) p->color_[0] << " " << (int) p->color_[1] << " " << (int) p->color_[2]
					<< " ";
			if (scale)
				pData << p->scale_3dx_ << " ";
			pData << std::endl;
		}
	}

	// now output the visibility
	if (visibility){
		for (const auto& p : patches) {
			if (binary) {
				uint32_t nrImgs = p->images_.size();
				pData.write((char*) &nrImgs, sizeof(uint32_t));
				for (uint32_t imgId : p->images_)
					pData.write((char*)&imgId, sizeof(uint32_t));
			} else {
				pData << (int) p->images_.size() << " ";
				for (int imgId : p->images_)
					pData << (uint32_t) imgId << " ";
				pData << std::endl;
			}
		}
	}
	pData.flush();
	pData.close();
}

template<typename Element>
void DynOctTree<Element>::toPly(const char* name, bool noScalar) {

	// collect the points and cells
	std::vector<Eigen::Vector3f> points;
	std::vector<Eigen::Vector3f> normals;
	std::vector<Eigen::Vector3i> colors;
	std::vector<float> pointlevel;
	std::vector<float> flatness;
	std::vector<Eigen::Vector4f> cells; // center and with
	std::vector<Eigen::Vector3i> cellColors;

	// now loop through the tree
	Leaf_iterator<Element> end = this->end();
	for (Leaf_iterator<Element> leaf = this->begin(); leaf != end; leaf++) {
		int cellLevel = nodeLevel(&(*leaf));
//		if (cellLevel != 8)
//			continue;
		std::vector<Element>& patches = leaf->data;
		Eigen::Vector3f cellColor(0, 0, 0);
		for (const auto& p : patches) {
			points.emplace_back(p->x(), p->y(), p->z());
			normals.emplace_back(p->normal_[0], p->normal_[1], p->normal_[2]);
			colors.emplace_back();

			int denom = 0;
			Eigen::Vector3f colorf = p->color_;
			colors.back()[0] = std::min(255, (int) floor(colorf[0] + 0.5f));
			colors.back()[1] = std::min(255, (int) floor(colorf[1] + 0.5f));
			colors.back()[2] = std::min(255, (int) floor(colorf[2] + 0.5f));

			cellColor += Eigen::Vector3f(colorf[0], colorf[1], colorf[2]);
			pointlevel.push_back(p->scale_3dx_);
			flatness.push_back(p->flatness_);
//			flatness.push_back(p->m_images.size()  );
		}

		if (patches.size() > 0) {
			cellColor /= patches.size();
			cellColors.emplace_back();
			cellColors.back()[0] = (int) std::min(255, (int) floor(cellColor[0] + 0.5f));
			cellColors.back()[1] = (int) std::min(255, (int) floor(cellColor[1] + 0.5f));
			cellColors.back()[2] = (int) std::min(255, (int) floor(cellColor[2] + 0.5f));
			cells.emplace_back(leaf->c_[0], leaf->c_[1], leaf->c_[2], leaf->width_);
		}

	}

	// output the points
	std::ofstream pFile((std::string(name) + "-points.ply").c_str(), std::ofstream::out);
	pFile << "ply" << '\n' << "format ascii 1.0" << '\n' << "element vertex " << (int) points.size()
			<< '\n' << "property float x" << '\n' << "property float y" << '\n'
			<< "property float z" << '\n' << "property float nx" << '\n' << "property float ny"
			<< '\n' << "property float nz" << '\n' << "property uchar red" << '\n'
			<< "property uchar green" << '\n' << "property uchar blue" << '\n';
	if (!noScalar)
		pFile << "property float scalar_scale" << '\n' << "property float scalar_flatness" << '\n';
	pFile << "end_header" << '\n';

	for (int ii = 0; ii < points.size(); ii++) {
		const auto& p = points[ii];
		const auto& c = colors[ii];
		const auto& n = normals[ii];
		pFile << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << n[0] << ' ' << n[1] << ' ' << n[2]
				<< ' ' << c[0] << ' ' << c[1] << ' ' << c[2];
		if (!noScalar)
			pFile << " " << pointlevel[ii] << " " << flatness[ii];
		pFile << '\n';
	}
	pFile.close();

	// output the tree visualization
	std::ofstream tFile(std::string(name) + "-tree.ply", std::ofstream::out);
	tFile << "ply" << '\n' << "format ascii 1.0" << '\n' << "element vertex "
			<< (int) 8 * cells.size() << '\n' << "property float x" << '\n' << "property float y"
			<< '\n' << "property float z" << '\n' << "property uchar diffuse_red" << '\n'
			<< "property uchar diffuse_green" << '\n' << "property uchar diffuse_blue" << '\n'
			<< "element face " << (int) 6 * cells.size() << '\n'
			<< "property list uchar int vertex_index" << '\n' << "end_header" << '\n';

	// the vertices
	int vertices[8][3] { 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0 };
	for (int ii = 0; ii < cells.size(); ii++) {
		float w = cells[ii][3];
		Eigen::Vector3f c(cells[ii][0] - w / 2, cells[ii][1] - w / 2, cells[ii][2] - w / 2);
		for (int v = 0; v < 8; v++) {
			tFile << c[0] + vertices[v][0] * w << " " << c[1] + vertices[v][1] * w << " "
					<< c[2] + vertices[v][2] * w << " " << cellColors[ii][0] << " "
					<< cellColors[ii][1] << " " << cellColors[ii][2] << "\n";
		}
	}

	// the faces
	for (int ii = 0; ii < cells.size(); ii++) {
		int si = 8 * ii;
		tFile << "4 " << si + 0 << " " << si + 1 << " " << si + 2 << " " << si + 3 << "\n";
		tFile << "4 " << si + 7 << " " << si + 6 << " " << si + 5 << " " << si + 4 << "\n";
		tFile << "4 " << si + 0 << " " << si + 4 << " " << si + 5 << " " << si + 1 << "\n";
		tFile << "4 " << si + 1 << " " << si + 5 << " " << si + 6 << " " << si + 2 << "\n";
		tFile << "4 " << si + 2 << " " << si + 6 << " " << si + 7 << " " << si + 3 << "\n";
		tFile << "4 " << si + 3 << " " << si + 7 << " " << si + 4 << " " << si + 0 << "\n";
	}

	tFile.flush();
	tFile.close();
}

template<typename Element>
void getBoundingBox(const std::vector<Element>& e, Eigen::Vector3f& min, Eigen::Vector3f& max) {
	// set unity cube if no patches
	if (e.size() == 0) {
		min = Eigen::Vector3f(-1, -1, -1);
		max = Eigen::Vector3f(1, 1, 1);
		return;
	}

	// evaluate min and max
	const float maxf = std::numeric_limits<float>::max();
	const float minf = std::numeric_limits<float>::min();

	min = Eigen::Vector3f(maxf, maxf, maxf);
	max = Eigen::Vector3f(minf, minf, minf);

	for (const auto& p : e) {
		min[0] = std::min(min[0], p->x());
		min[1] = std::min(min[1], p->y());
		min[2] = std::min(min[2], p->z());

		max[0] = std::max(max[0], p->x());
		max[1] = std::max(max[1], p->y());
		max[2] = std::max(max[2], p->z());
	}
}

#endif
