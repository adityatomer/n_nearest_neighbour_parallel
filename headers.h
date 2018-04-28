#ifndef HEADERS_H
#define HEADERS_H

#include<iostream>
#include<cstdio>
#include<cmath>
#include<list>
#include<vector>
#define DIMENSIONS 3
#define BOXEDGES 8
#define _sqr(x) ((x)*(x))
using namespace std;

class TreeNode;

struct PointCord{
	double cord[DIMENSIONS];
	TreeNode *node;

	public:
		PointCord(){

		}

		PointCord(double x, double y, double z){
			cord[0]=x;
			cord[1]=y;
			cord[2]=z;
		}
};


// struct point{
// 	double cord[DIMENSIONS];
// 	TreeNode *node;

// 	Point(){

// 	}

// 	Point(double x, double y, double z){
// 			cord[0]=x;
// 			cord[1]=y;
// 			cord[2]=z;
// 	}
// };


struct Box{
	PointCord boxEdges[BOXEDGES];
	PointCord center;
	double radius;
	Box(){}
	Box(double x, double xdim, double y, double ydim, double z, double zdim) {
		int k =0;
		boxEdges[k++] = PointCord(x, y, z);
		boxEdges[k++] = PointCord(x+xdim, y, z);
		boxEdges[k++] = PointCord(x+xdim, y+ydim, z);
		boxEdges[k++] = PointCord(x, y+ydim, z);    
		boxEdges[k++] = PointCord(x, y, z+zdim);
		boxEdges[k++] = PointCord(x+xdim, y, z+zdim);
		boxEdges[k++] = PointCord(x+xdim, y+ydim, z+zdim);
		boxEdges[k++] = PointCord(x, y+ydim, z+zdim);
		radius = sqrt(_sqr(xdim*0.5)+_sqr(ydim*0.5)+_sqr(zdim*0.5));
		center = PointCord(x+xdim/2, y+ydim/2, z+zdim/2);
  	}

  	Box(PointCord &p){
  		center=PointCord(p.cord[0],p.cord[1],p.cord[2]);
  		radius=0;
  		for(int i=0;i<BOXEDGES;++i){
  			boxEdges[i]=center;
  		}
  	}

};


struct PointSetPart{
	int start;
	int length;
	PointSetPart(int start, int length){
		this->start=start;
		this->length=length;
	}

};

class TreeNode{
	public:
		TreeNode* left;
		TreeNode* right;
		TreeNode* parent;
		Box boudingbox;
		int point_idx;
		std::list<TreeNode*> *interacts;

		TreeNode(TreeNode* parent){
			this->parent=parent;
			this->interacts=NULL;
		}
		TreeNode(TreeNode* parent,int point_idx){
			this->parent=parent;
			this->point_idx=point_idx;
			this->left=NULL;
			this->right=NULL;
			this->interacts=NULL;
		}

		TreeNode(TreeNode* left, TreeNode* right, TreeNode* parent, int point_idx){
			this->left=left;
			this->right=right;
			this->parent=parent;
			this->point_idx=point_idx;
			this->interacts=NULL;
		}

		void addInteracts(TreeNode *node){
			if(this->interacts==NULL){
				this->interacts=new std::list<TreeNode*>();
			}
			this->interacts->push_back(node);
		}
};

#endif