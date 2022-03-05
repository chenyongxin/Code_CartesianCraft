// Controller.cc
// Controller for the programme.
// Code_CartesianCraft (c) Yongxin Chen

#include "Controller.h"

Controller::Controller():
  _probe_index(MAX_PROBES, 3),
  _probe_mask (MAX_PROBES   ),
  _slice_index(MAX_SLICES, 3),
  _slice_mask (MAX_SLICES, 3)
{}

void Controller::init(const MyMPI &mympi){

  // Correct BC
  for(int d=0; d<3; d++){
    if(!mympi.lboundary(d)){
      _left_vbc[d] = NONE;
      _left_sbc[d] = NONE;
      _left_pbc[d] = NONE;
    }

    if(!mympi.rboundary(d)){
      _right_vbc[d] = NONE;
      _right_sbc[d] = NONE;
      _right_pbc[d] = NONE;
    }
  }

  // Init timer and make initial recording
  t0 = std::chrono::steady_clock::now(); 
  t_last = t0;
  tick(mympi);
}

void Controller::set_probe(const int i, const Vec3d &coords, const Grid &grid, const MyMPI &mympi){
  if(i >= MAX_PROBES){
    if(mympi.rank() == 0) {cerr << "Probe index exceeds the capacity." << endl;}
    exit(1);
  }
  if( grid.inside(coords, mympi) ){
    _probe_mask(i) = true;
    Vec3i index = grid.coords2index(coords);
    for(int d=0; d<3; d++) {_probe_index(i, d) = index(d);}
  }
}

void Controller::set_slice(const int i, const int d, const double coord, const Grid &grid, const MyMPI &mympi){
  if(i >= MAX_SLICES){
    if(mympi.rank() == 0) {cerr << "Slice index exceeds the capacity." << endl;}
    exit(1);
  }
  if( (grid.lrange(d) <= coord) && (coord < grid.rrange(d)) ){
    _slice_mask(i, d) = true;
    for(int j=grid.start(d); j<grid.end(d); j++){
      if( (grid.x(d, j) <= coord) && (coord < grid.x(d, j+1)) ){
        _slice_index(i, d) = j;
        break;
      }
    }
  }
}

void Controller::tick(const MyMPI &mympi){

  // Get time difference and update.
  t_now = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds  = t_now-t0;
  std::chrono::duration<double> timestep_seconds = t_now-t_last;
  t_last = t_now;

  // Output result 
  if(mympi.rank() == 0){
    int width = 25;
    ofstream file;
    file.open(MONITOR_TIMER, ios::app);
    if( _step % INTERVALS == 0 ){
      file << setw(width) << "Steps";
      file << setw(width) << "Time_step";
      file << setw(width) << "Simulation_time";
      file << setw(width) << "Time_per_step_(s)";
      file << setw(width) << "Elapsed_time_(s)";
      file << setw(width) << "Elapsed_time_(h:m:s)";
      file << endl;
    }
    double t = elapsed_seconds.count();
    int    h = t/3600;
    int    m = t/60-h*60;
    int    s = t-h*3600-m*60;
    string time_breakdown = to_string(h).append(":").append(to_string(m).append(":").append(to_string(s)));
    file << setw(width) << _step;
    file << setw(width) << _dt;
    file << setw(width) << _time;
    file << setw(width) << timestep_seconds.count();
    file << setw(width) << t;
    file << setw(width) << time_breakdown;
    file << endl; 
    file.close();
  }
}
