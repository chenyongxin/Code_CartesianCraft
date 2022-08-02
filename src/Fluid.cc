// Fluid.cc
// Top class of the solver.
// Code_CartesianCraft (c) Yongxin Chen

#include "Fluid.h"
#ifdef USE_HDF
#include <mpi.h>
#include "hdf5.h"
#endif

using namespace std;

void Fluid::init(){

  // Init objects and make directories
  initialiser.read_file(mympi, controller, grid);
  initialiser.init(mympi, controller, grid);
  psolver.init(mympi, grid, controller);
  bc.init(mympi, controller, grid);
  ibm.init(grid, controller, mympi);
  stat.init(grid);

  // Init fields
  u.allocate (grid.n());
  h.allocate (grid.n()); 
  p.allocate (grid.n());
  ed.allocate(grid.n());
  r.allocate(grid.n());
  for(int d=0; d<3; d++){
    uij[d].allocate(grid.n());
    uij_cen[d].allocate(grid.n());
  }
  div.allocate(grid.n());
  cfl.allocate(grid.n());
  if(controller.get_dragmodel() != 0){
    D.allocate(grid.n());
  }

  // Assign variables
  nu = controller.nu();

  // Read fields if restart
  if(initialiser.restart()){
    
    // Read fields and controls
    read();
  
    // Resume time-averaged statistics 
    if(initialiser.restatistics()){
#ifdef USE_HDF
      stat.read_hdf(mympi, grid);
#else
      stat.read(mympi); 
#endif
    }
  }
  else{
    // new run, apply initial condition.
    bc.apply_vic(u, controller, grid);
  } // end restart

  bc.apply_vbc(u, mympi, controller, grid);

  // Initial step velocity gradient tensor and eddy viscosity 
  Discretisation::gradtensor(uij, uij_cen, u, grid);
  SGS::compute_SGS(ed, u, uij_cen, grid, controller);
  bc.apply_sbc(ed, mympi, controller);
}

void Fluid::updateAB(){

  // Turbulent inlet 
  if(controller.turbulent_inflow_generation()){
    if((controller.step() % TURBINLET_INTERVALS) == 0) {bc.fluctuation_gen(mympi, controller, grid);}
  }

  Discretisation::conv_diff(r, u, grid, nu, ed, uij);
  
  // Get drag force term
  if(controller.get_dragmodel() != 0){
    Discretisation::get_dragforce(u, D, grid, controller.Cd(), controller.cdheight());
  }
  
  if(controller.step() == 0) {h = r;} // It will reduce the init step to Euler
  Discretisation::AB(u, r, h, controller.bforce(), controller.dt());
  
  // Add additional terms
  if(controller.get_dragmodel() != 0){
    Discretisation::add_vectorfields(u, D);
  }

  h = r;                              // Update history term
  updateUP(u);
  
  // Get new velocity gradient tensor and eddy viscosity based on the new velocitiy
  Discretisation::gradtensor(uij, uij_cen, u, grid);
  SGS::compute_SGS(ed, u, uij_cen, grid, controller);
  bc.apply_sbc(ed, mympi, controller);
}

void Fluid::updateUP(VectorFieldd &v){
  bc.apply_vbc(v, mympi, controller, grid);
  if(controller.get_ibmtype() != 0) {
    ibm.ibc(v);
    bc.apply_vbc(v, mympi, controller, grid);
  }
  Discretisation::divergence(div, v, grid);
  div.scale(-1./controller.dt());
  psolver.solve(div, p, mympi);
  bc.apply_pbc(p, mympi, controller);
  Discretisation::project(v, p, grid, controller);
  bc.apply_vbc(v, mympi, controller, grid);
}

double Fluid::cfl_limit(){
  cfl = 0.0;    // reset du/dx
  for(int d=0; d<3; d++){
    for(int k=u.ks(); k<u.ke(); k++){
      for(int j=u.js(); j<u.je(); j++){
        for(int i=u.is(); i<u.ie(); i++){
          Vec3i ijk(i,j,k);
          cfl(i,j,k) += abs(u[d](i,j,k) * grid.dxi(d,ijk(d)));
        }
      }
    }
  }
  return  controller.cfl() / cfl.max(mympi); // get the minimum time step
}

void Fluid::flows(){
  init();
  
  // Iterate 
  while(controller.time() < controller.totaltime()){

    // Update an adaptive time step
    if(controller.adaptive()) controller.set_dt( cfl_limit() );

    // Time integration
    updateAB();

    // Increment timer
    controller.increment_time();
    controller.increment_step();
    controller.tick(mympi);
    if(mympi.rank() == 0) {cout << controller.step() << endl;}

    // Whole field visualisation
    field_visu();

    // Slice visualisation
    slice_visu();

    // Do statistics   
    time_average();

    // IO for restart
    if( (controller.step() > 1) && (controller.step()%50000 == 0) ){ save(); } 
  }
}

void Fluid::read(){
 
#ifdef USE_HDF
  // Collectively read fields
  hid_t file_id = H5Fopen( (string(HDF_DIR)+string("/")+string("Fluid.h5")).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  H5read_field(file_id, "u0", mympi, grid, u[0]);
  H5read_field(file_id, "u1", mympi, grid, u[1]);
  H5read_field(file_id, "u2", mympi, grid, u[2]);
  H5read_field(file_id, "h0", mympi, grid, h[0]);
  H5read_field(file_id, "h1", mympi, grid, h[1]);
  H5read_field(file_id, "h2", mympi, grid, h[2]);
  H5read_field(file_id, "p",  mympi, grid, p   );

  H5Fclose(file_id);

  // Read controls
  ifstream file;
  file.open( (string(SAVE_DIR)+string("/")+string("Fluid")).c_str(), ios::binary );

  int i;
  double d; 
  file.read((char*)&i, sizeof(int));         // read time step
  file.read((char*)&d, sizeof(double));      // read simulation time
  controller.set_step(i);                    
  controller.set_time(d);

  file.read((char*)&d, sizeof(double));      // read Delta_t  
  controller.set_dt(d);                      

  file.read((char*)&i, sizeof(int));         // read visu step
  file.read((char*)&d, sizeof(double));      // read visu time
  visu_step = i;
  controller.set_visu_time(d); 

  file.read((char*)&i, sizeof(int));         // read slice step
  file.read((char*)&d, sizeof(double));      // read slice time
  slice_step = i;   
  controller.set_slice_time(d); 
   
  file.close();
#else
  // Open a file 
  ifstream file;
  char name[100];
  sprintf(name, "%s/Fluid_x%dx%dx%d", SAVE_DIR, mympi.coords(0), mympi.coords(1), mympi.coords(2));  
  file.open(name, ios::binary);
  
  // Read controls
  int i;
  double d; 
  file.read((char*)&i, sizeof(int));         // read time step
  file.read((char*)&d, sizeof(double));      // read simulation time
  controller.set_step(i);                    
  controller.set_time(d);
  
  file.read((char*)&d, sizeof(double));      // read Delta_t  
  controller.set_dt(d);                      

  file.read((char*)&i, sizeof(int));         // read visu step
  file.read((char*)&d, sizeof(double));      // read visu time
  visu_step = i;
  controller.set_visu_time(d); 

  file.read((char*)&i, sizeof(int));         // read slice step
  file.read((char*)&d, sizeof(double));      // read slice time
  slice_step = i;   
  controller.set_slice_time(d); 
   
  // Read fields
  u.read(file);
  h.read(file);
  p.read(file);

  // Close a file 
  file.close();
#endif
}

void Fluid::save() const {

#ifdef USE_HDF
  // Collectively write fields
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate( (string(HDF_DIR)+string("/")+string("Fluid.h5")).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  H5Pclose(fapl);

  H5write_field(file_id, "u0", mympi, grid, u[0]);
  H5write_field(file_id, "u1", mympi, grid, u[1]);
  H5write_field(file_id, "u2", mympi, grid, u[2]);
  H5write_field(file_id, "h0", mympi, grid, h[0]);
  H5write_field(file_id, "h1", mympi, grid, h[1]);
  H5write_field(file_id, "h2", mympi, grid, h[2]);
  H5write_field(file_id, "p",  mympi, grid, p   );

  H5Fclose(file_id);

  // Write controls
  if(mympi.master()){
    ofstream file;
    file.open( (string(SAVE_DIR)+string("/")+string("Fluid")).c_str(), ios::binary );
    
    int i = controller.step();
    double d = controller.time();
    file.write((char*)&i, sizeof(int));         // write time step
    file.write((char*)&d, sizeof(double));      // write simulation time

    d = controller.dt();
    file.write((char*)&d, sizeof(double));      // write Delta_t 

    i = visu_step;
    d = controller.visu_time();
    file.write((char*)&i, sizeof(int));         // write visu step
    file.write((char*)&d, sizeof(double));      // write visu time

    i = slice_step;
    d = controller.slice_time();
    file.write((char*)&i, sizeof(int));         // write slice step
    file.write((char*)&d, sizeof(double));      // write slice time

    file.close();
  }

  // Save statistics
  stat.save_hdf(mympi, grid);
#else
  // Open a file 
  ofstream file;
  char name[100];
  sprintf(name, "%s/Fluid_x%dx%dx%d", SAVE_DIR, mympi.coords(0), mympi.coords(1), mympi.coords(2));  
  file.open(name, ios::binary);
  
  // Write controls
  int i = controller.step();
  double d = controller.time();
  file.write((char*)&i, sizeof(int));         // write time step
  file.write((char*)&d, sizeof(double));      // write simulation time
  
  d = controller.dt();
  file.write((char*)&d, sizeof(double));      // write Delta_t

  i = visu_step;
  d = controller.visu_time();
  file.write((char*)&i, sizeof(int));         // write visu step
  file.write((char*)&d, sizeof(double));      // write visu time

  i = slice_step;
  d = controller.slice_time();
  file.write((char*)&i, sizeof(int));         // write slice step
  file.write((char*)&d, sizeof(double));      // write slice time

  // Write fields
  u.save(file);
  h.save(file);
  p.save(file);

  // Close a file
  file.close();

  // Save Statistics. Only save the data. Whether to read statistics depends on controller.
  stat.save(mympi);
#endif
}

void Fluid::field_visu(){
  if(controller.time() >= controller.visu_time()){
    int total_digits = 4;
    int n_reserve    = 5;
    visu_step++;
    vector<Fieldd> additional_fields;
    vector<string> additional_fieldnames;
    additional_fields.reserve(n_reserve);
    additional_fieldnames.reserve(n_reserve);

    // Add additional fields and their names
    additional_fields.push_back(ed);
    additional_fieldnames.push_back(string("nu_sgs"));

#ifdef USE_HDF
    H5fluid_fields(u, p, additional_fields, additional_fieldnames, grid, mympi, (string(HDF_DIR)+string("/")).c_str(),
        string(total_digits-to_string(visu_step).length(), '0').append(to_string(visu_step)).c_str());
#else
    myvtk.fluid_fields(u, p, additional_fields, additional_fieldnames, grid, mympi, (string(VTK_DIR)+string("/")).c_str(),
        string(total_digits-to_string(visu_step).length(), '0').append(to_string(visu_step)).c_str());
#endif
    controller.increment_visu();

    // Write file name and its simulation time to file.
    if(mympi.rank() == 0){
      ofstream file;
      file.open(MONITOR_VISU, ios::app);
      file << string(total_digits-to_string(visu_step).length(), '0').append(to_string(visu_step)) << "   ";
      file << controller.time() << endl;
      file.close();
    }
  } // end if current time greater than visualisation time
}

void Fluid::slice_visu() {
  if(controller.time() >= controller.slice_time()){
    int total_digits = 4;
    slice_step++;
    void *dummy;
    myvtk.slice_fields(u, p, (Fieldd *) dummy, (char **) &dummy, 0, grid, mympi, controller, (string(SLICE_DIR)+string("/")).c_str(),
        string(total_digits-to_string(visu_step).length(), '0').append(to_string(visu_step)).c_str());
    controller.increment_slice();

    // Write file name and its simulation time to file.
    if(mympi.rank() == 0){
      ofstream file;
      file.open(MONITOR_SLICE, ios::app);
      file << string(total_digits-to_string(slice_step).length(), '0').append(to_string(slice_step)) << "   ";
      file << controller.time() << endl;
      file.close();
    }
  }
}

void Fluid::time_average(){
  if(controller.time() >= controller.stat_time()) {
    stat.do_statistics(u, p); 
  }
  if(controller.step()%100000 == 0) {stat_visu();}
}

void Fluid::stat_visu(){
  stat.write(grid, mympi);
}

Fluid::~Fluid(){
  stat_visu();
  save();
}
