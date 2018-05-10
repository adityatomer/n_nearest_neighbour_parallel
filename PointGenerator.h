#ifndef _DEF_POINTGEN_
#define _DEF_POINTGEN_

#include <vector>
#include <math.h>

struct double3Coord {
    double coord[3];
    double& operator[] (unsigned int idx) {
            return coord[idx];
    }
    double3Coord () {}
    double3Coord (double x, double y, double z) {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }
};

         
double
normalize(double3Coord& vector)
{
    double dist = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    if (dist > 1e-6)
    {
        vector[0] /= dist;
        vector[1] /= dist;
        vector[2] /= dist;
    }
    return dist;
}

float 
dot(double3Coord v0, double3Coord v1)
{
    return v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2];
}

double3Coord
cross(double3Coord v0, double3Coord v1)
{
    double3Coord rt;
    rt[0] = v0[1]*v1[2]-v0[2]*v1[1];
    rt[1] = v0[2]*v1[0]-v0[0]*v1[2];
    rt[2] = v0[0]*v1[1]-v0[1]*v1[0];	
    return rt;
}


void createShellPoints(double clusterScale, int numBodies, std::vector<double3Coord> & PointCoords) {
    // Ripped from CUDA NBody demo
    double scale = clusterScale;
    double inner = 2.5 * scale;
    double outer = 4.0 * scale;
    
    int i = 0;
    while (i < numBodies)//for(int i=0; i < numBodies; i++) 
    {
        double x, y, z;
        x = rand() / (double) RAND_MAX * 2 - 1;
        y = rand() / (double) RAND_MAX * 2 - 1;
        z = rand() / (double) RAND_MAX * 2 - 1;
        
        double3Coord ptt    = double3Coord(x, y, z);
        double len = normalize(ptt);
        if (len > 1)
            continue;
        
        double3Coord pt;
        pt[0] =    ptt[0] * (inner + (outer - inner) * rand() / (double) RAND_MAX);
        pt[1] =    ptt[1] * (inner + (outer - inner) * rand() / (double) RAND_MAX);
        pt[2] =    ptt[2] * (inner + (outer - inner) * rand() / (double) RAND_MAX);
        
        PointCoords.push_back(pt);
        i++;
    }
}

void createRandomPoints(double clusterScale, double velocityScale, int numBodies,    std::vector<double3Coord> & PointCoords) {
    // Note ripped from CUDA n-body demo!
    double scale = clusterScale * std::max(1.0, numBodies / (1024.0));
    double vscale = velocityScale * scale;
    
    int i = 0;
    while (i < numBodies) 
    {
                double3Coord pt;
		//const int scale = 16;
		pt[0] = rand() / (double) RAND_MAX * 2 - 1;
		pt[1] = rand() / (double) RAND_MAX * 2 - 1;
		pt[2] = rand() / (double) RAND_MAX * 2 - 1;        
     
		pt[0] *= scale; // pos[0]
                pt[1] *= scale; // pos[0]
                pt[2] *= scale; // pos[0]
        
	        i++; 
                PointCoords.push_back(pt);
    }
}

void createGridPoints(int width, int height, std::vector<double3Coord> & PointCoords) {
    for(int x=0; x<width; x++) {
        for(int y=0; y<height; y++) {
            double3Coord pt;
            pt[0] = x;
            pt[1] = y;
            pt[2] = 0.0;    
            PointCoords.push_back(pt);
        }
    }
}
 

#endif
