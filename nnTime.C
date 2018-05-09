#include <iostream>
#include <algorithm>
#include "gettime.h" 
// #include "utils.h"
// #include "cilk.h" 
#include <stdlib.h>
#include <sys/time.h>
#include <iomanip>
#include <iostream>
#include "parseCommandLine.h"
#include "PointGenerator.h"
#include "NNeighbour.h"
#include <vector>
#include "NNeighbour.h"
#include<assert.h>
#define CILK 0
#define SQR(x) ((x)*(x))
using namespace std;

NNeighbour * nneighbourSet = NULL;
int _nPointCoords;

double distsqr(double3 a, double3 b) {
    return (SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]));
}

int nearest_neighbor_naive(std::vector<double3> & PointCoords, int orig) {
    double minsqr=1e30; 
    int minidx=-1;
    for(int i=0; i<PointCoords.size(); i++) {
        if (i != orig) {
            double r2 = distsqr(PointCoords[i], PointCoords[orig]);
            if (r2<minsqr) {
                minidx = i;
                minsqr = r2;
            }
        }
    }
    assert(minidx>=0 && minidx != orig);
    return minidx;
}

void addPoint(int i, double3 p) {
    nneighbourSet->addPoint(i, p[0],p[1],p[2]);
}

// Used for correctness check.
double3 getNearestPoint(int PointCoordidx) {
    PointCoord p = nneighbourSet->getPoints()[nneighbourSet->getPoints()[PointCoordidx].nn];
    return double3(p.coord[0], p.coord[1], p.coord[2]);
}

bool checkNN(vector<double3> & PointCoordset) {
    // Use max 0.5 secs for checking or max 1000 PointCoords
    int maxPointCoords = 1000;
    double maxtime = 0.5;
    timer ctimer;
    ctimer.start();
    srand(time(NULL));
    for(int i=0; i<maxPointCoords; i++) {
        if (ctimer.total() < maxtime || i < 50) {
            int PointCoordcheck = random()%PointCoordset.size();
            double3 nearestPointCoord = PointCoordset[nearest_neighbor_naive(PointCoordset, PointCoordcheck)];
            double3 alg_nn = getNearestPoint(PointCoordcheck);
            // Two PointCoords may be same distance, so it is not enough to check correct PointCoord index
            if (distsqr(PointCoordset[PointCoordcheck], alg_nn) == distsqr(PointCoordset[PointCoordcheck], nearestPointCoord)) {
                // ok
            } else {
                std::cout.precision(5);
                std::cout << "Fail, algorithm returned incorrect nn for PointCoord " << PointCoordcheck << std::endl;
                std::cout << "Distance to alg's nn: " << distsqr(PointCoordset[PointCoordcheck], alg_nn) << " correct: " <<  
                        distsqr(PointCoordset[PointCoordcheck], nearestPointCoord) << std::endl;
                std::cout << "Failed query PointCoord was: " << PointCoordset[PointCoordcheck][0] << " " << PointCoordset[PointCoordcheck][1] << " "   
                        << PointCoordset[PointCoordcheck][2] << std::endl;
                std::cout << "Algorithm returned: " << alg_nn[0] << " " << alg_nn[1] << " " << alg_nn[2] << std::endl;
                return false;
            }
        } else {
            std::cout << ":: Exceeded 0.5 sec checking limit." << std::endl;
            break;
        }
    }
    return true;
}

void computeNN() {
    nneighbourSet->recompute();
    nneighbourSet->compute_wsr_nn();
}

void createNNeighbourInstance(size_t numOfPoints) {
    if ( nneighbourSet != NULL ) delete(nneighbourSet);
    nneighbourSet = new NNeighbour();
    nneighbourSet->reserve(numOfPoints);
    _nPointCoords = numOfPoints;
}
 
void computeNearestNeighbour(std::string testType, size_t numPoints, char * outFile, int rounds) {
    vector<double3> PointCoordset;
  
    srand(1);
    if (testType == "grid") {
        createGridPoints((int) sqrt(numPoints), (int) sqrt(numPoints), PointCoordset);
    }
    else if (testType == "random") {
        createRandomPoints(10.0, 10.0, numPoints, PointCoordset);
    }
    else if (testType == "shell") {
        createShellPoints(10.0,  numPoints, PointCoordset);
    } else assert(false);

    for (int i=0; i < rounds; i++) {
        startTime();
        createNNeighbourInstance(PointCoordset.size());
        for(int i=0; i<PointCoordset.size(); i++) {
            addPoint(i, PointCoordset[i]);
        }

        computeNN();
        nextTimeN();
        
        bool success = checkNN(PointCoordset);
        if (outFile != NULL) {
            FILE * outf = fopen(outFile, "w");
            fprintf(outf, "%d\n", (success ? 1 : 0));
            fclose(outf);
        }
     }
     cout << endl;    
}

int main(int argc, char* argv[]) {
    commandLine P(argc,argv," [-o <outFile>] [-r <rounds>] [-n <numPoints>] [-t <testType=random|shell|grid>]");
    char* testType = P.getOptionValue("-t");
    char* oFile = P.getOptionValue("-o");
    int rounds = P.getOptionIntValue("-r",1);
    size_t numPoints = (size_t) P.getOptionIntValue("-n", 1000);
    computeNearestNeighbour(std::string(testType), numPoints, oFile, rounds);
}
