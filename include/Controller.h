// Controller.h
// Controller for the programme.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Vec.h"
#include "Grid.h"
#include "Array.h"
#include "MyMPI.h"
#include "parameter.h"

using std::ios;
using std::setw;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::to_string;

class Controller{

  private:

    // Timing
    int    _step        =    0;  // time step
    double _dt          = -1.0;  // time step Delta t
    double _time        =  0.0;  // simulation time
    double _walltime    = -1.0;  // wall time
    double _runtime     =  0.0;  // elapsed run time
    double _totaltime   = -1.0;  // total simulation time   
    bool   _adaptive    = true;  // adaptive time stepping
    double _cfl         =  0.5;  // CFL limit

    std::chrono::time_point<std::chrono::steady_clock> t0;     // t0 of simulation
    std::chrono::time_point<std::chrono::steady_clock> t_last; // time in the last time step
    std::chrono::time_point<std::chrono::steady_clock> t_now;  // time in this time step

    // Fluid properties
    double _nu   = -1.;          // kinematic viscosity
    double _utau = -1.;          // u_\tau
    double _z0   = -1.;          // z0
    double _uref =  1.;          // reference velocity
    double _Cd   =  0.;          // fluid drag coefficient, resulting from the unresolved canopy element
    double _cdheight = 0.;       // canopy height, under where the drag is applied 
      
    // Cases
    int _case    = -1;           // case for ic and bc types 
    int _ibmtype =  0;           // IBM roughness type
    int _dragmodel = 0;          // Use drag model for the unresolved elements. 0 for off and 1 for on

    // Body force
    Vec3d _bforce;

    // Boundary conditions
    // Dummy velocity
    int  _left_vbc[3] = {INLET,     PERIODIC, SYMMETRY};
    int _right_vbc[3] = {OUTFLOW,   PERIODIC, SYMMETRY};

    // Pressure
    int  _left_pbc[3] = {DIRICHLET, PERIODIC, NEUMANN };
    int _right_pbc[3] = {NEUMANN,   PERIODIC, NEUMANN };

    // Scalar
    int  _left_sbc[3] = {DIRICHLET, PERIODIC, NEUMANN };
    int _right_sbc[3] = {NEUMANN,   PERIODIC, NEUMANN };

    // Turbulent inflow generation
    bool _turbulent_inflow_generation = false;

    // Probes
    Array1b _probe_mask;
    Array2i _probe_index;

    // Slices
    Array2b _slice_mask;
    Array2i _slice_index;

    // Visualisation for the whole field
    double _visu_dt   = 1.e8;       // Time step increment for the next visualisation event
    double _visu_time = 1.e8;       // Time to do visualisation. Can be used as t0 and timer

    // Visualisation for slice field
    double _slice_dt   = 1.e8;      // Time step increment for the next visualisation event
    double _slice_time = 1.e8;      // Time to do visualisation. Can be used as t0 and timer

    // Time average statistics 
    double _stat_time  = 1.e8;      // Time average statistics initialisation time

    // Initial velocity condition
    vector<Vec3d> _vic_coords_start; // velocity initial condition coordinate start
    vector<Vec3d> _vic_coords_end;   // velocity initial condition coordinate end
    vector<Vec3d> _vic_value;        // velocity initial condition values

  public:

    // Contructor and deconstructor
    Controller();
    ~Controller(){}

    // Set boundary conditions
    void set_left_vbc (const int d, const int x){  _left_vbc[d] = x; }
    void set_right_vbc(const int d, const int x){ _right_vbc[d] = x; }
    void set_left_pbc (const int d, const int x){  _left_pbc[d] = x; }
    void set_right_pbc(const int d, const int x){ _right_pbc[d] = x; }
    void set_left_sbc (const int d, const int x){  _left_sbc[d] = x; }
    void set_right_sbc(const int d, const int x){ _right_sbc[d] = x; }

    // Get boundary conditions
    int left_vbc (const int d) const { return  _left_vbc[d]; }
    int right_vbc(const int d) const { return _right_vbc[d]; }
    int left_pbc (const int d) const { return  _left_pbc[d]; }
    int right_pbc(const int d) const { return _right_pbc[d]; }
    int left_sbc (const int d) const { return  _left_sbc[d]; }
    int right_sbc(const int d) const { return _right_sbc[d]; }

    // Initalise controller. Called after BC are set
    // by correcting MPI interface BC types to NONE and init timer.
    void init(const MyMPI &);

    // Set timing
    void set_step     (const int    x) {_step     = x;}
    void set_dt       (const double x) {_dt       = x;}
    void set_time     (const double x) {_time     = x;}
    void set_walltime (const double x) {_walltime = x;}
    void set_totaltime(const double x) {_totaltime= x;}
    void set_adaptive (const bool   x) {_adaptive = x;}
    void set_cfl      (const double x) {_cfl      = x;}

    // Get timing
    int    step()      const {return _step;     }
    double dt()        const {return _dt;       }
    double time()      const {return _time;     }
    double walltime()  const {return _walltime; }
    double totaltime() const {return _totaltime;}
    bool   adaptive()  const {return _adaptive; }
    double cfl()       const {return _cfl;      }
    double runtime()   const {return _runtime;  }

    // Reference velocity
    void   set_uref    (const double x) { _uref = x;  }
    double uref        ()         const {return _uref;}

    // Visualisation
    // Set values
    void set_visu_dt   (const double x) {_visu_dt   = x;}
    void set_visu_time (const double x) {_visu_time = x;}
    void set_slice_dt  (const double x) {_slice_dt  = x;}
    void set_slice_time(const double x) {_slice_time= x;}
    // Get values
    double visu_time()  const {return _visu_time; }
    double slice_time() const {return _slice_time;}
    // Increment for the next event
    void increment_visu () {_visu_time  += _visu_dt; }
    void increment_slice() {_slice_time += _slice_dt;}

    // Time average statistics
    void   set_stat_time(const double x) {_stat_time = x;}
    double stat_time()   const {return _stat_time;}

    // Set body force
    void set_bforce   (const int d, const double x) {_bforce(d) = x;}

    // Get body force
    double bforce     (const int d) const {return _bforce(d);}
    Vec3d  bforce     ()            const {return _bforce;   }

    // Set IBM roughness type
    void set_ibmtype (const int x) {_ibmtype = x;}

    // Get IBM roughness type
    int get_ibmtype() const {return _ibmtype; }

    // Set drag model for unresolved elements, switch on
    void set_dragmodel(const int x) {_dragmodel = x;}

    // Get drag model status
    int get_dragmodel() const {return _dragmodel;}

    // Increment timing
    void increment_time () {   _time  += _dt;} // increment simulation time by adding time step
    void increment_step () {   _step  +=   1;} // increment simulation step by adding 1
    void tick(const MyMPI &);                  // record iteration timing

    // Set ith probe with a specific coordinate
    void set_probe(const int i, const Vec3d &coords, const Grid &grid, const MyMPI &mympi);

    // Set ith slice location with a specific coordinate 
    void set_slice(const int i, const int d, const double coord, const Grid &grid, const MyMPI &mympi);

    // Check if probe and slice is inside the partition. Loop with parameter MAX_..
    bool probe_inside(const int i)               const {return _probe_mask(i)   ;}
    bool slice_inside(const int i, const int d)  const {return _slice_mask(i, d);}

    // Get probe and slice index
    Vec3i probe_index(const int i)               const {return Vec3i(_probe_index(i,0),
                                                                     _probe_index(i,1),
                                                                     _probe_index(i,2));}

    int slice_index (const int i, const int d)   const {return _slice_index(i, d);}

    // Set and get case type
    void set_case   (const int    x)                   {_case    = x;}
    int  get_case   ()                           const {return _case;}

    // Set and get fluid properties
    void   set_nu   (const double x)                   {_nu      = x;}
    double nu       ()                           const {return _nu;  }

    void   set_utau (const double x)                   {_utau    = x;}
    double utau     ()                           const {return _utau;}

    void   set_z0   (const double x)                   {_z0      = x;}
    double z0       ()                           const {return   _z0;}

    void   set_Cd   (const double x)                   {_Cd      = x;}
    double Cd       ()                           const {return _Cd;  }

    void   set_cdheight(const double x)                {_cdheight = x;}
    double cdheight ()                           const {return _cdheight;}
    
    // Set and get turbulence inflow generation.
    void set_turbulent_inflow_generation(const bool x) {_turbulent_inflow_generation = x;}
    bool turbulent_inflow_generation() const {return _turbulent_inflow_generation;}

    // Set and get velocity initial condition.
    void set_vic(const Vec3d &start, const Vec3d &end, const Vec3d &value){
      _vic_coords_start.push_back(start);
      _vic_coords_end.push_back(end);
      _vic_value.push_back(value);
    }
    Vec3d vic_coords_start(const int i) const {return _vic_coords_start[i];}
    Vec3d vic_coords_end  (const int i) const {return _vic_coords_end[i];  }
    Vec3d vic_value       (const int i) const {return _vic_value[i];       }
    int   vic_number      ()            const {return _vic_value.size();   }
};

#endif // CONTROLLER_H
