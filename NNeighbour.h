#ifndef __NNEIGHBOUR_
#define __NNEIGHBOUR_

#include <cstdio>
#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <assert.h>
#include <cmath>
#define DIMENSIONS 3
#define BOXEDGES 8
#define MAXPARALLELDEPTH 5
#define _sqr(x) ((x)*(x))


class TreeNode;
struct PointCoord {
    double coord[DIMENSIONS];
    double m;
    int nn;
    TreeNode * node;

    PointCoord(){ 
        nn = (-1); 
    } 
    PointCoord(double x, double y, double z) {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;    
        nn = (-1); 
    }
    double rsqr(PointCoord &p) {
        assert(DIMENSIONS == 3);
        return _sqr(p.coord[0]-coord[0])+_sqr(p.coord[1]-coord[1])+_sqr(p.coord[2]-coord[2]);
    } 
};

struct Box {
    PointCoord edges[BOXEDGES];
    PointCoord center;
    double r;
    Box() {

    }
    Box(double x, double xdim, double y, double ydim, double z, double zdim) {
        int k =0;
        edges[k++] = PointCoord(x, y, z);
        edges[k++] = PointCoord(x+xdim, y, z);
        edges[k++] = PointCoord(x+xdim, y+ydim, z);
        edges[k++] = PointCoord(x, y+ydim, z);    
        edges[k++] = PointCoord(x, y, z+zdim);
        edges[k++] = PointCoord(x+xdim, y, z+zdim);
        edges[k++] = PointCoord(x+xdim, y+ydim, z+zdim);
        edges[k++] = PointCoord(x, y+ydim, z+zdim);
        r = sqrt(_sqr(xdim*0.5)+_sqr(ydim*0.5)+_sqr(zdim*0.5));
        center = PointCoord(x+xdim/2, y+ydim/2, z+zdim/2);
    }
  
    Box(PointCoord & pt) {
        r = 0;
        center = PointCoord(pt.coord[0], pt.coord[1], pt.coord[2]);
        for(int k=0; k<BOXEDGES; k++) edges[k] = center;
    }
};

class TreeNode {
    public:
        TreeNode * left;
        TreeNode * right;
        TreeNode * parent;
        Box boundingBox;
        int depth;
        double lengthMax;
        double totalmass;
        PointCoord masscenter;
        std::list<TreeNode *> * interacts;
        int PointCoord_idx;

        TreeNode(TreeNode * p) { 
            parent = p; interacts = NULL;  
        }
        TreeNode(TreeNode *  l, TreeNode *  r, TreeNode *  p) { 
            left = l; right = r; parent = p;interacts = NULL;   
        }
        TreeNode(TreeNode *  l, TreeNode *  r, TreeNode *  p, int  pt_idx) { 
            left = l; right = r; parent = p; 
            PointCoord_idx = pt_idx;  interacts = NULL; 
        }

        ~TreeNode() {
            if (interacts != NULL) delete(interacts);
        }

        void add_interaction(TreeNode * nd) {
            if (interacts == NULL){
                interacts = new std::list<TreeNode *>();
            };
            interacts->push_back(nd);
        }
};

struct PointCoordset_part {
    int begin;
    int length;
    PointCoordset_part(int b, int l) { begin = b; length = l; } 
};

class NNeighbour {
    public:
    
    NNeighbour(): s(2.01) {
        PointCoordidx = NULL;
        workarray = NULL;
        internal_nodes = NULL;
        leaves = NULL;
        root = NULL;
    }
  
    ~NNeighbour() {
        free(PointCoordidx);
        free(workarray);
        free(leaves);
        free(internal_nodes);
    }
  
    void reserve(int n) {
        PointCoords.resize(n);
    } 
  
    void addPoint(int i, double x, double y, double z) {
        PointCoords[i] = PointCoord(x,y,z);
    }
   
    // Return the min and max value of a PointCoord for given dimension
    std::pair<double, double> PointCoordset_dim(int dim) {
        std::vector<PointCoord>::iterator iter;
        double dmin = 1e30, dmax=-1e30;
        for(iter=PointCoords.begin(); iter!=PointCoords.end(); ++iter) {
            dmin = std::min(dmin, iter->coord[dim]);
            dmax = std::max(dmax, iter->coord[dim]);
        }
        return std::pair<double, double>(dmin, dmax);
    }
  
    std::vector<PointCoord> & getPoints() {
        return PointCoords;
    }
  
    void recompute() {
        Boxes.clear();
        buildKDTree();
        buildWSR();
    }
  
    std::vector<Box> & getBoxes() {
        return Boxes;
    }
  
    std::vector< std::pair<PointCoord, PointCoord> > &  getPairs() {
        return pairs;
    }
  
  
    void buildWSR() {
        pairs.clear();
        wellSeparatedRegion(root);
    }
  
  
    int nearest_neighbor_naive(int pntidx) {
        PointCoord & p = PointCoords[pntidx];
        double minDistSqr = 1e30;
        int minDistIdx = -1;
        for(int i=PointCoords.size()-1; i>=0; i--) {
            if (i != pntidx) {
                double r2 = p.rsqr(PointCoords[i]);
                if (r2 < minDistSqr) {
                    minDistSqr = r2;
                    minDistIdx = i;
                }
            }
        }
        return minDistIdx;
    }

    void inorder(TreeNode *root){
        if(root==NULL){
            return;
        }

        inorder(root->left);
        cout<<root<<" "<<root->PointCoord_idx<< "\n";
        inorder(root->right);
    }

    void preorder(TreeNode *root){
        if(root==NULL){
            return;
        }
        cout<<root<<" "<<root->PointCoord_idx<< "\n";
        inorder(root->left);
        inorder(root->right);
    }
  
    void printPoint(PointCoord p){
        cout<<"("<<p.coord[0]<<", "<<p.coord[1]<<", "<<p.coord[2]<<")\n";
    }

    void compute_nn_naive() {
        int pointCoordSize = PointCoords.size();
        #ifdef CILK
        cilk_for(int i=0; i<pointCoordSize; i++) {
        #else
        for(int i=0; i<pointCoordSize; i++) {
        #endif
            PointCoords[i].nn = nearest_neighbor_naive(i); 
        }
    }
    
    std::pair<int,double> mindist(PointCoord & p, TreeNode * node, int minDistIdx, double minDistSqr){ 
        if (node->left == NULL) {
            double r2 = p.rsqr(PointCoords[node->PointCoord_idx]);
            if (r2<minDistSqr) 
                return std::pair<int,double>(node->PointCoord_idx, r2);
            else 
                return std::pair<int,double>(minDistIdx, minDistSqr);
        } else {
            std::pair<int,double> minleft = mindist(p, node->left, minDistIdx, minDistSqr);
            return mindist(p, node->right, minleft.first, minleft.second);
        }
    }
    
    // Check only via interacting edges
    int nearest_neighbor(int pntidx) {
        return PointCoords[pntidx].nn;
    }
    
    void compute_wellSeparatedRegion_nn() {
        int pointCoordSize = PointCoords.size();
        #ifdef CILK
            cilk_for(int pntidx=0; pntidx<pointCoordSize; pntidx++){ 
        #else 
            for(int pntidx=0; pntidx<pointCoordSize; pntidx++){
        #endif
                double minDistSqr = 1e30;
                int minDistIdx = -1;
                PointCoord & p = PointCoords[pntidx];
                std::list<TreeNode *> & interacts = *(p.node->interacts);
                std::list<TreeNode *>::iterator iter;
                //std::cout << "Interacts: " << pntidx << " " << interacts.size() << std::endl;
                for(iter = interacts.begin(); iter != interacts.end(); iter++) {
                    std::pair<int, double> mind = mindist(p, *iter, minDistIdx, minDistSqr);
                    minDistIdx = mind.first;
                    minDistSqr = mind.second;
                    //std::cout << pntidx << " : " <<    minDistIdx << " : " << minDistSqr << std::endl;
                }
                assert(minDistIdx >= 0);
                p.nn = minDistIdx;
        }
    }

        
private:
    std::vector<PointCoord> PointCoords;
    std::vector< std::pair<PointCoord, PointCoord> > pairs;
    
    int * PointCoordidx;
    int * workarray;
    std::vector<Box> Boxes;
    TreeNode * internal_nodes;
    TreeNode * leaves;
    TreeNode * root;
    double s;


    void printPointNode(){
        for(int i=0;i<PointCoords.size();++i){
            // cout<<"i: "<<i<<" "<<" Node: "<<PointCoords[i].node<<" => ";
            printPoint(PointCoords[i]);
            cout<<"\n";
        }
    }
    
    // Building the tree (using the simpler method, so it is easier to parallize)
    void buildKDTree() {
        int tsz = 0;
        PointCoordidx = (int*) realloc(PointCoordidx, PointCoords.size()*sizeof(int));
        workarray = (int *) realloc(workarray, PointCoords.size()*sizeof(int));
        leaves = (TreeNode *) realloc(leaves, PointCoords.size()*sizeof(TreeNode));
        internal_nodes = (TreeNode *) realloc(internal_nodes, PointCoords.size()*sizeof(TreeNode));
        for(int i=0; i<PointCoords.size(); i++) {
                PointCoordidx[i] = i;
        }
        root = kdTree(PointCoordset_part(0, PointCoords.size()), NULL, 0);
    }
    
    // Nodes are preallocated to avoid doing mallocs (as new()).
    // Leaves have their own array, leaves are indexed by their PointCoord index.
    TreeNode *    leaf(TreeNode *    parent, int ptidx) {
        TreeNode * node = &leaves[ptidx];
        *node = TreeNode(NULL, NULL, parent, ptidx);
        PointCoords[ptidx].node = node;
        node->boundingBox = Box(PointCoords[ptidx]);
        node->lengthMax = 0;
        if (parent == NULL) {
             node->depth = 0;
        } else {
             node->depth = parent->depth+1;
        }
        return node;
    }
    
    // Internal nodes are indexed by the pivot position when a PointCoordset
    // is split.
    TreeNode * newnode(TreeNode * parent, int splitpos) {
        TreeNode * node = &internal_nodes[splitpos];
        *node = TreeNode(parent);
                
         if (parent == NULL) {
            node->depth = 0;
         } else {
            node->depth = parent->depth+1;
         }
            
        return node;
    }
    
    std::pair<double, double> PointCoordpart_dim(PointCoordset_part part, int dim) {
        double dmin = 1e30, dmax=-1e30;
        for(int i=0; i<part.length; i++) {
            PointCoord & pt = PointCoords[PointCoordidx[i+part.begin]];
            dmin = std::min(dmin, pt.coord[dim]);
            dmax = std::max(dmax, pt.coord[dim]);
        }
        return std::pair<double, double>(dmin, dmax);
    }
    
    Box PointCoordpart_Box(PointCoordset_part part) {
        std::pair<double,double> xd = PointCoordpart_dim(part,0);
        std::pair<double,double> yd = PointCoordpart_dim(part,1);
        std::pair<double,double> zd = PointCoordpart_dim(part,2);
        return Box(xd.first,xd.second-xd.first, yd.first, yd.second-yd.first, zd.first, zd.second-zd.first); 
    }
    
    double lengthMax(Box bb) {
        double mx = 0;
        for(int i=1; i<BOXEDGES; i++) {
            for(int d=0; d<DIMENSIONS; d++) {
                mx = std::max(mx, fabs(bb.edges[i].coord[d]-bb.edges[i-1].coord[d]));
            }
        }
        return mx;
    }
    
    std::pair<int, double> maxpart_dim(PointCoordset_part part) {
        double maxdim = 0;
        double maxdimmin = 0;
        int maxdimnum = -1;
        for(int i=0; i<DIMENSIONS; i++) {
            std::pair<double,double> dm = PointCoordpart_dim(part, i);
            if (dm.second-dm.first >= maxdim) {
                maxdim = dm.second-dm.first;
                maxdimmin = dm.first;
                maxdimnum = i;
            }
        }
        return std::pair<int,double>(maxdimnum, maxdimmin + maxdim*0.5);
    }
    
    TreeNode * kdTree(PointCoordset_part p, TreeNode * parent, int level) {
        assert(p.length>=1);
        if (p.length == 1) {
            return leaf(parent, PointCoordidx[p.begin]);
        } else {
            std::pair<int,double> dmax = maxpart_dim(p);
            // Split
            double midPointCoord = dmax.second;
            int end = p.begin + p.length;
            int c=p.begin;
            for(int i=p.begin; i<end; i++) {
                if (PointCoords[PointCoordidx[i]].coord[dmax.first] < midPointCoord) {
                    workarray[c++] = PointCoordidx[i];
                }
            }
            int p2start=c;
            for(int i=p.begin; i<end; i++) {
                if (PointCoords[PointCoordidx[i]].coord[dmax.first] >= midPointCoord) {
                    workarray[c++] = PointCoordidx[i];
                }
            }
            // Move to correct parts of the array
            for(int i=p.begin; i<end; i++) {
                PointCoordidx[i] = workarray[i];
            }
            TreeNode * node = newnode(parent, p2start);
            node->boundingBox = PointCoordpart_Box(p);
            node->lengthMax = lengthMax(node->boundingBox);
    
            #ifdef CILK
                TreeNode * leftn;
                TreeNode * rightn;
                if (level < MAXPARALLELDEPTH) {
                        leftn = cilk_spawn kdTree(PointCoordset_part(p.begin, p2start-p.begin), node, level+1);
                        rightn =     kdTree(PointCoordset_part(p2start, c-p2start), node, level+1);
                        cilk_sync;
                } else {
                        rightn = kdTree(PointCoordset_part(p2start, c-p2start), node, level+1);
                        leftn = kdTree(PointCoordset_part(p.begin, p2start-p.begin), node, level+1);
                }
                node->left = leftn;
                node->right = rightn;
            #else
                node->left = kdTree(PointCoordset_part(p.begin, p2start-p.begin), node, level+1);
                node->right = kdTree(PointCoordset_part(p2start, c-p2start), node, level+1);
            #endif
            // For visualization:
            //Boxes.push_back(tree[newnodeid].boundingBox);
            return node;
        }
    }
    
    bool wellSeparated(Box & b1, Box & b2) {
        // Smallest radius of a ball that can contain aither rectangle
        double r = std::max(b1.r, b2.r);
        double d = sqrt(b1.center.rsqr(b2.center))-2*r;
        return (d>s*r);
    }
    
    bool wellSeparated(PointCoord & p, Box & b2, double s) {
        double r = b2.r;
        double d = sqrt(b2.center.rsqr(p))-2*r;
        return (d>s*r);
        
    }
    
    
    void wellSeparatedRegion(TreeNode * node) {
        if (node->left == NULL) {
            // is leaf
            return;
        } else {
            #ifdef CILK
                cilk_spawn wellSeparatedRegion(node->left);
                wellSeparatedRegion(node->right);
                cilk_sync;
            #else
                wellSeparatedRegion(node->left);
                wellSeparatedRegion(node->right);
            #endif
            wellSeparatedRegionPair(node->left, node->right);
        }
    }
    
    void wellSeparatedRegionPair(TreeNode * t1, TreeNode * t2) {
        if (wellSeparated(t1->boundingBox, t2->boundingBox)) {
            if (t1->left == NULL) t1->add_interaction(t2);
            if (t2->left == NULL) t2->add_interaction(t1);
            //pairs.push_back(std::pair<PointCoord,PointCoord> (t1.boundingBox.center, t2.boundingBox.center));
        } else if ((t1->lengthMax) > (t2->lengthMax)) {
            wellSeparatedRegionPair(t1->left, t2);
            wellSeparatedRegionPair(t1->right, t2);
        } else {
            wellSeparatedRegionPair(t1, t2->left);
            wellSeparatedRegionPair(t1, t2->right);
        }
    }
    
};

#endif
