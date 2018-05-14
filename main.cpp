#include <iostream>
#include <algorithm>
#include "gettime.h" 
// #include "utils.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <sys/time.h>
#include <unistd.h>
#include<chrono>
using namespace std;
#include "PointGenerator.h"
#include "NNeighbour.h"
#include <ctype.h>
#include <stdio.h>
#include <vector>
#include<fstream>
#include<assert.h>
#define SQUARE(x) ((x)*(x))


NNeighbour * nneighbourSet = NULL;
int _nPointCoords;

double distsqr(double3Coord a, double3Coord b) {
    return (SQUARE(a[0]-b[0])+SQUARE(a[1]-b[1])+SQUARE(a[2]-b[2]));
}

int nearest_neighbor_naive(std::vector<double3Coord> & PointCoords, int orig) {
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

void addPoint(int i, double3Coord p) {
    nneighbourSet->addPoint(i, p[0],p[1],p[2]);
}

// Used for correctness check.
double3Coord getNearestPoint(int PointCoordidx) {
    PointCoord p = nneighbourSet->getPoints()[nneighbourSet->getPoints()[PointCoordidx].nn];
    return double3Coord(p.coord[0], p.coord[1], p.coord[2]);
}

bool checkNN(vector<double3Coord> & PointCoordset) {
    // Use max 0.5 secs for checking or max 1000 PointCoords
    int maxPointCoords = 1000;
    double maxtime = 100000;
    timer ctimer;
    ctimer.start();
    srand(time(NULL));
    for(int i=0; i<maxPointCoords; i++) {
        if (ctimer.total() < maxtime || i < 50) {
            int PointCoordcheck = random()%PointCoordset.size();
            double3Coord nearestPointCoord = PointCoordset[nearest_neighbor_naive(PointCoordset, PointCoordcheck)];
            double3Coord alg_nn = getNearestPoint(PointCoordcheck);
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
    nneighbourSet->compute_wellSeparatedRegion_nn();
}

void createNNeighbourInstance(size_t numOfPoints) {
    if ( nneighbourSet != NULL ) delete(nneighbourSet);
    nneighbourSet = new NNeighbour();
    nneighbourSet->reserve(numOfPoints);
    _nPointCoords = numOfPoints;
}


void writeToFile(std::string path,vector<double3Coord> & PointCoordset){
    ofstream myfile(path);
    for(int i=0;i<PointCoordset.size();++i){
        int PointCoordcheck=i;
        double3Coord p                = PointCoordset[i];
        PointCoord p_neighbourPC = nneighbourSet->getPoints()[nneighbourSet->getPoints()[PointCoordcheck].nn];
        double3Coord p_neighbour= double3Coord(p_neighbourPC.coord[0], p_neighbourPC.coord[1], p_neighbourPC.coord[2]);
        if (myfile.is_open()){
            myfile <<"("<<p.coord[0]<<", "<<p.coord[1]<<", "<<p.coord[2]<<")";
            myfile <<"  =>  ";
            myfile <<"("<<p_neighbour.coord[0]<<", "<<p_neighbour.coord[1]<<", "<<p_neighbour.coord[2]<<")";
            myfile <<"  =>  "<<distsqr(p, p_neighbour);
            myfile <<"\n";
        }
    }
    myfile.close();
}
 
void computeNearestNeighbour(std::string testType, size_t numPoints, std::string outFile, int rounds) {
    vector<double3Coord> PointCoordset;
  
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
        writeToFile(outFile,PointCoordset);
        bool success = checkNN(PointCoordset);
        if(!success){
            cout<<"ERROR!!!";
        }
     }
     cout << endl;    
}

void commandLineArgs(int argc,char *argv[], std::string &oFile, int &rounds, size_t &numPoints,std::string &testType, int &workers){
    int c=0;
    while ((c = getopt (argc, argv, "w::o:r::n::t:")) != -1)
    switch (c)
    {   
        case 'o':
            oFile=optarg;
            break;
        case 'r':
            rounds=atoi(optarg);
            break;
        case 'n':
            numPoints=atoi(optarg);
            break;
        case 't':
            testType = optarg;
            break;
        case 'w':
            workers = atoi(optarg);
            break;
        case '?':
            if (optopt == 'c')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
        default:
            return;
    }
}

int main(int argc, char* argv[]) {
    std::string testType;
    std::string oFile;
    int rounds = 1;
    size_t numPoints = 10;
    int workers=1;
    commandLineArgs(argc,argv,oFile, rounds, numPoints, testType,workers);
    #ifdef CILK
       __cilkrts_set_param("nworkers", argv[2]);  
    #endif
    __int64 now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();	
    

    computeNearestNeighbour("random", atoi(argv[1]), "./output.out", 1);    
    
    __int64 now1 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    float timeTaken=float (now1- now);
    std::cout<<"time taken: "<<timeTaken/1000<<"\n";

    //computeNearestNeighbour(std::string(testType), numPoints, oFile, rounds);
    return 1;
}
