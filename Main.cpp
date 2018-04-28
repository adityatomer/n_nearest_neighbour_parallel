#include<iostream>
#include "headers.h"



// void createGridPoints(int width, int height, std::vector<double3> & points) {
//   for(int x=0; x<width; x++) {
//     for(int y=0; y<height; y++) {
//       double3 pt;
//       pt[0] = x;
//       pt[1] = y;
//       pt[2] = 0.0;  
//       points.push_back(pt);
//     }
//   }
// }


std::pair<double,double>min_max_point_dim(vector<PointCord>points, PointSetPart p, int dim){
	std::pair<double,double>res;
	double maxi=0;
	double mini=1e30;
	for(int i=0;i < points.size();++i){
		maxi=max(maxi, points[i].cord[dim]);
		mini=min(maxi, points[i].cord[dim]);
	}

	return std::pair<double,double>(mini,maxi);
}

std::pair<int, double> maxpart_dim(vector<PointCord>points, PointSetPart part) {
	int maxDimNum=0;
	double maxDimVal=0;
	double maxDimStart=0;

	for(int i=0;i<DIMENSIONS;++i){
		std::pair<double,double>minmax= min_max_point_dim(points,i);

		if(minmax.second-minmax.first > maxDimVal){
			maxDimVal=minmax.second - minmax.first;
			maxDimNum=i;
			maxDimStart=minmax.first;
		}
	}
	return std::pair<double,double>(maxDimNum,maxDimStart+maxDimVal*0.5);
}


Box boundingBox(vector<PointCord>points,PointSetPart pp){
	
	std::pair<double,double>x=min_max_point_dim(points,pp,0);
	std::pair<double,double>y=min_max_point_dim(points,pp,1);
	std::pair<double,double>z=min_max_point_dim(points,pp,2);

	Box b=Box(x.first,x.second-x.first, y.first,y.second-y.first, z.first,z.second-z.first);
	return b;
}

TreeNode* leaf(vector<PointCord>points, TreeNode *parent, int pointIndx){
	TreeNode tn=new TreeNode(parent,pointIndx);
	points[pointIndx].node=tn;

	node->boudingbox=box(points[pointIndx]);
	if(parent==NULL){
		node->level=0;
	}else{
		node->level=parent->level+1;
	}
	return node;
}


TreeNode* newnode(TreeNode *parent, ){
	TreeNode* node=new TreeNode(parent);
	if(parent==NULL){
		node->level=0;
	}else{
		node->level=parent->level+1;	
	}
	return node;
}


//treenode * btree(pointset_part p, treenode * parent, int level) {
TreeNode* buildTreeInternal(vector<PointCord>points, pointset_part p,treenode *parent, int level){
	//change this
	assert(p.length>=1);
	if (p.length == 1) {
      return leaf(points, parent, p.begin);
    }
    std::pair<int,double>dmax=maxpart_dim(points, p);
    // 
    int mid=dmax.second;
	int en=p.start+p.length;
	int c=0;
    for(int i=0;i<points.size();++i){
    	if(points[i] < mid){
    		c++;
    	}
    }
    TreeNode *tn=newnode(parent,c);
    tn->boudingbox=boundingBox(points,p);

    tn->left=buildTreeInternal(vector<PointCord>points, pointset_part(p.begin,c-p.begin),tn,level+1);
    tn->right=buildTreeInternal(vector<PointCord>points, pointset_part(mid,en-c),tn,level+1);

    return tn;
}

void inorder(TreeNode *root){
	if(root==NULL){
		return;
	}

	inorder(root->left);
	cout<<" inorder => "<<node->level;
	inorder(root->right);
}

void buildTree(vector<PointCord>points) {
	pointset_part p=pointset_part(0, points.size());
	TreeNode *root = buildTreeInternal(points, p, NULL, 0);


}


int main(){

	// NNeighbour nn=NNeighbour();
	vector<PointCord>pointset;
	
	int testType=0;//GRID

	if(testType==0){
		int gridSize=4;
		int length=gridSize;
		int width=gridSize;
		int height=0;

		for(int i=0;i<length;++i){
			for(int j=0;j<width;++j){
				PointCord pc=PointCord(i,j,0.0);
				pointset.push_back(pc);
			}

		}

		for(int i=0;i<gridSize*gridSize;++i){
			cout<<pointset[i]<<" ";
		}
		cout<<endl<<endl;
	}



	buildTree(pointset);


	return 0;
}