#incldue<iostream>
#include<cstdio>
#include<cmath>
#define DIMENSIONS 3
#define BOXEDGES 8
#define _sqr(x) ((x)*(x))
using namespace std;

struct PointCord{
	double cord[DIMENSIONS];

	public:
		PointCord(){

		}

		PointCord(double x, double y, double z){
			cord[0]=x;
			cord[1]=y;
			cord[2]=z;
		}
};

struct point{
	double cord[DIMENSIONS];

	Point(){

	}

	Point(double x, double y, double z){
			cord[0]=x;
			cord[1]=y;
			cord[2]=z;
	}
};

struct Box{
	PointCord[BOXEDGES] boxEdges;
	PointCord center;
	double radius;
	box(){

	}
	box(double x, double xdim, double y, double ydim, double z, double zdim) {
		int k =0;
		edges[k++] = staticpoint(x, y, z);
		edges[k++] = staticpoint(x+xdim, y, z);
		edges[k++] = staticpoint(x+xdim, y+ydim, z);
		edges[k++] = staticpoint(x, y+ydim, z);    
		edges[k++] = staticpoint(x, y, z+zdim);
		edges[k++] = staticpoint(x+xdim, y, z+zdim);
		edges[k++] = staticpoint(x+xdim, y+ydim, z+zdim);
		edges[k++] = staticpoint(x, y+ydim, z+zdim);
		r = sqrt(_sqr(xdim*0.5)+_sqr(ydim*0.5)+_sqr(zdim*0.5));
		center = PointCord(x+xdim/2, y+ydim/2, z+zdim/2);
  	}
};

class TreeNode{

	public:
		TreeNode* left;
		TreeNode* right;
		TreeNode* parent;
		box boudingbox;
		int point_idx;
		std::list<TreeNode*> *interacts;


};