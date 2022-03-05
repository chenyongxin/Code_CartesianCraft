// Initialiser.h
// Initialise control parameters and flow fields
// Code_CartesianCraft (c) Yongxin Chen

#ifndef INITIALISER_H
#define INITIALISER_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

#include "BC.h"
#include "Vec.h"
#include "Grid.h"
#include "MyMPI.h"
#include "parameter.h"
#include "Controller.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;
using std::istringstream;

class Initialiser{

  private:

    // Temporal parameters saved for init objects.
    // Grid
    int grid_types[3] = {0,0,0};

    // MPI
    int mpi_blocks[3] = {0,0,0};

    // New run or re-start
    bool _restart = false;

    // New time-averaged statistics or resume it
    bool _restatistics = false;
      
    // Probe
    vector<int>    probe_index;
    vector<Vec3d>  probe_coords;

    // Slice
    vector<int>    slice_index;
    vector<int>    slice_direction;
    vector<double> slice_coords;

  public:

    Initialiser(){}
    ~Initialiser(){}

    // Initailise controller and MPI by reading parameters from input file.
    void read_file(const MyMPI &, Controller &, Grid &);

    // Initialise controller, MPI and grid and IBM.
    void init(MyMPI &, Controller &, Grid &) const;

    // Enquiry if resume running the case.
    bool restart() const {return _restart;}
  
    // Enquiry if resume time-averaged statistics
    bool restatistics() const {return _restatistics;}
};
#endif // INITIALISER_H
