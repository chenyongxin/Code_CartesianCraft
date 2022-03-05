// parameter.h
// Definition of parameters in the programme.
// Code_CartesianCraft (c) Yongxin Chen
// DO NOT CHANGE unless you know how to do so.

#ifndef PARAMETER_H
#define PARAMETER_H

// Ghost cell
#define GC 2

// Constant
#define PI          3.141592653589793
#define TWO_PI      6.283185307179586
#define C_VREMAN    0.1          // for Cs=0.2 
//#define C_VREMAN  0.07225      // for Cs=0.17
#define aND         0.2          // for stability
#define BETA        0.1          // Gamma scheme coefficient
#define MAGNITUDE   0.05         // initial white noise magnitude

// Boundary condition
#define NONE        0
#define PERIODIC    1
#define NEUMANN     2
#define DIRICHLET   3
#define INLET       4
#define OUTFLOW     5
#define CONVECTIVE  6
#define SYMMETRY    7
#define WALL        8        // no-slip wall
#define TURBINLET   9        // turbulent boundary inlet bc
#define WALLMODEL   10       // wall function

// Cases for initial and boundary conditions
//#define CASE_USER      0
#define CASE_TBL         1
#define CASE_URBAN       2
#define CASE_CYLINDER    3
#define CASE_OPENCHANNEL 4

// IBM roughness type
#define CUBOIDS 1            // homogeneous cuboids
#define URBAN   2            // urban roughness

// Folder name
#define   HDF_DIR   "HDF"    // save HDF data
#define  GRID_DIR   "GRID"   // save grid data
#define  DATA_DIR   "DATA"   // save all data 
#define  SAVE_DIR   "SAVE"   // save checkpoint data
#define  STAT_DIR   "STAT"   // save time averaged statistics VTK data
#define SLICE_DIR   "SLICE"  // save slice VTK data
#define   VTK_DIR   "VTK"    // save whole flow field VTK data

// File name
#define MONITOR_HYPRE "MONITOR_Psolver"
#define MONITOR_TIMER "MONITOR_Timer"
#define MONITOR_VISU  "MONITOR_Visu"
#define MONITOR_SLICE "MONITOR_Slice"

// Maximum passive scalars
#define MAX_SCALARS 2

// Maximum slices for visualisation in each direction
#define MAX_SLICES  10

// Maximum probes in fluid
#define MAX_PROBES  50

// Timing record intervals
#define INTERVALS 10

// Turbulent inflow generation intervals, an integer number greater than 1
// can reduce the time in fluctuation generation.
#define TURBINLET_INTERVALS 5

// HYPRE solver
#define HYPRE_TOLERANCE     1.0e-6
#define HYPRE_AMG_PRINT_LV  0
#define HYPRE_PCG_PRINT_LV  0
#define HYPRE_MAX_ITERATION 100

#endif // PARAMETER_H
