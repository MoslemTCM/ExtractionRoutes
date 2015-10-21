/**
 * Filtre Bresenam realise dans le cadre du cours SYS809 de ETS
 * Auteurs: Matthieu CAUFFIEZ, Benjamin MILLUY
 *
 */


#ifndef PYRAMIDTREE_H_
#define PYRAMIDTREE_H_


#include <list>
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"

using namespace std;

typedef unsigned char CharPixelType; // IO
const unsigned int Dimension = 2;
typedef itk::Image<CharPixelType, Dimension> ImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

struct beamletData
{
	ImageType::IndexType start;
	ImageType::IndexType end;
	int lenth;
};

struct FortChildsTree
{
	ImageType::RegionType::IndexType regionStart;
	ImageType::RegionType::SizeType regionSize;
	list<beamletData> beamletlist;
	FortChildsTree *upperLeft;
	FortChildsTree *upperRight;
	FortChildsTree *lowerLeft;
	FortChildsTree *lowerRight;
	FortChildsTree *parent;

struct FortChildsTree* newFortChildsTree(ImageType::RegionType::IndexType regionStart,
ImageType::RegionType::SizeType regionSize)
{
	FortChildsTree *node = new FortChildsTree;
	node->regionStart = regionStart;
	node->regionSize = regionSize;
	node->beamletlist.clear();
	node->upperLeft = NULL;
	node->upperRight = NULL;
	node->lowerLeft = NULL;
	node->lowerRight = NULL;
	node->parent = NULL;
	return node;
}

struct FortChildsTree* insertTreeNode(FortChildsTree *node, ImageType::RegionType::IndexType regionStart,
		ImageType::RegionType::SizeType regionSize, short int position)
{
	FortChildsTree *newNode = new FortChildsTree;
	newNode->regionStart = regionStart;
	newNode->regionSize = regionSize;
	newNode->beamletlist.clear();
	newNode->upperLeft = NULL;
	newNode->upperRight = NULL;
	newNode->lowerLeft = NULL;
	newNode->lowerRight = NULL;
	newNode->parent = node;

	switch (position)
	{
	case 1: {
		node->upperLeft = newNode;
		return newNode;
	}
	case 2:{
		node-> upperRight = newNode;
		return newNode;
	}
	case 3:{
		node->lowerLeft = newNode;
		return newNode;
	}
	case 4:{
		node->lowerRight = newNode;
		return newNode;
	}
	}

}

void  getLevelElements(struct FortChildsTree *root, int level, std::list<FortChildsTree*> *listOfTrees  ,int iteration=0 )
{
if(level == 0)
{
	return;
}

if(iteration != level )
{
	getLevelElements(root->upperLeft, level, listOfTrees, iteration + 1);
	getLevelElements(root->upperRight, level, listOfTrees ,iteration + 1);
	getLevelElements(root->lowerLeft, level, listOfTrees, iteration + 1);
	getLevelElements(root->lowerRight, level, listOfTrees, iteration + 1);
}else
{
	listOfTrees->push_back(root);
}


}
};
#endif /* PYRAMIDTREE_H_ */
