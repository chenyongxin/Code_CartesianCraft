// Fluid.h
// Top class of the solver.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef FLUID_H
#define FLUID_H

#include <cmath>
#include <vector>
#include <string>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "BC.h"
#include "Vec.h"
#include "IBM.h"
#include "SGS.h"
#include "Grid.h"
#include "Stat.h"
#include "H5IO.h"
#include "MyVTK.h"
#include "MyMPI.h"
#include "Field.h"
#include "Psolver.h"
#include "parameter.h"
#include "Controller.h"
#include "Initialiser.h"
#include "VectorField.h"
#include "Discretisation.h"

using std::ios;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::to_string;
using std::vector;

class Fluid{

  private:

    // Necessary objects
    BC           bc;
    IBM          ibm;
    Grid         grid;
    MyMPI        mympi;
    MyVTK        myvtk;
    Psolver      psolver;
    Controller   controller;
    Initialiser  initialiser;
    Stat         stat;

    // Primary fields
    VectorFieldd u;                 // velocities
    VectorFieldd h;                 // history terms for Adam-Bashforth
    Fieldd       p;                 // pressure
    Fieldd       ed;                // eddy viscosity

    // Assistant fields	          
    Fieldd       div;               // divergence of velocity
    VectorFieldd r;                 // conv-diff term
    VectorFieldd uij[3];            // staggered velocity gradient tensor
    VectorFieldd uij_cen[3];        // velocity gradient tensor in cell-centre
    Fieldd       cfl;               // compute local cfl number
    VectorFieldd D;                 // drag force due to unresolved element

    // Variables		          
    double       nu;                // kinematic viscosity
    int          visu_step =0;      // whole field visualisation step
    int          slice_step=0;      // slice visualisation step

    // Functions		          
    void init    ();                // initialisation
    void updateAB();                // Adam-Bashforth time integration
    void updateUP(VectorFieldd&);   // update velocity U and pressure P
    void read    ();                // read all fluid and statistics data
    void save    () const;          // save all fluid and statistics data for restart
    double cfl_limit();             // compute time step based on the CFL limit
    void field_visu();              // visualise the whole field
    void slice_visu();              // visualise a slice
    void time_average();            // do time average statitics 
    void stat_visu();               // visualise statistics result

  public:
    Fluid(int &argc, char **argv): mympi(argc, argv) {}
    ~Fluid();

    void flows();
};

#endif // FLUID_H
