// Initialiser.cc
// Initialise control parameters and flow fields
// Code_CartesianCraft (c) Yongxin Chen

#include "Initialiser.h"

// Make string character lower case
void lower_string(string &str){
  for(int i=0; i<str.length(); i++){
    if(str[i] >= 'A' && str[i] <= 'Z') { str[i] += 32; }
  }
}

void Initialiser::read_file(const MyMPI &mympi, Controller &controller, Grid &grid) {
  ifstream file("input.ccc");
  if(!file.is_open()){
    if(mympi.rank()==0) { cerr << "Input file ``input.ccc'' does not exist." << endl; }
    exit(1);
  }
  string line;
  while(getline(file, line)){
    istringstream str(line);
    string word, keyword, list;
    str >> word;
    if(word[0] == '$' && word[word.length()-1] == '$') {
      keyword = word.substr(1, word.length()-2);
      lower_string(keyword);    // get command name
      getline(file, line);      // get command list
      istringstream list(line); // wrap to input stream 
      if(keyword == "mpi"){
        for(int d=0; d<3; d++) list >> mpi_blocks[d];
      }
      else if(keyword == "grid"){
        int d, type;
        list >> d >> type;
        grid_types[d] = type;
        if(type == 1){
          double neg, pos, h;
          list >> neg >> pos >> h;
          grid.set_uniform_line(d, neg, pos, h);
        }
        else if(type == 2){
          double uneg, upos, pos, h, factor;
          list >> uneg >> upos >> pos >> h >> factor;
          grid.set_stretch_line(d, uneg, upos, pos, h, factor);
        }
        else{
          if(mympi.rank()==0) cerr << "Unknown type when initialising grid.";
          exit(1);
        }
      }
      else if(keyword == "case"){
        int icase;
        list >> icase;
        controller.set_case(icase);
      }
      else if(keyword == "utau"){
        double utau;
        list >> utau;
        controller.set_utau(utau);
      }
      else if(keyword == "z0"){
        double z0;
        list >> z0;
        controller.set_z0(z0);
      }
      else if(keyword == "nu"){
        double nu;
        list >> nu;
        controller.set_nu(nu);
      }
      else if(keyword == "inflow"){
        int x;
        list >> x;
        controller.set_turbulent_inflow_generation(x);
      }
      else if(keyword == "totaltime"){
        double totaltime;
        list >> totaltime;
        controller.set_totaltime(totaltime);
      }
      else if(keyword == "walltime"){
        double walltime;
        list >> walltime;
        controller.set_walltime(walltime);
      }
      else if(keyword == "dt"){
        double dt;
        list >> dt;
        controller.set_dt(dt);
      }
      else if(keyword == "adaptive"){
        int adaptive;
        list >> adaptive;
        controller.set_adaptive(adaptive);
      }
      else if(keyword == "cfl"){
        double cfl;
        list >> cfl;
        controller.set_cfl(cfl);
      }
      else if(keyword == "bodyforce"){
        int d;
        double x;
        list >> d >> x;
        controller.set_bforce(d, x);
      } 
      else if(keyword == "probe"){
        int i;
        double x, y, z;
        list >> i >> x >> y >> z;
        probe_index.push_back(i);
        probe_coords.push_back(Vec3d(x,y,z));
      }
      else if(keyword == "slice"){
        int i, d;
        double coord;
        list >> i >> d >> coord;
        slice_index.push_back(i);
        slice_direction.push_back(d);
        slice_coords.push_back(coord);
      }
      else if(keyword == "restart"){
        list >> _restart;	
      }
      else if(keyword == "restatistics"){
        list >> _restatistics;
      }
      else if(keyword == "ibmtype"){
        int type;
        list >> type;
        controller.set_ibmtype(type);
      }
      else if(keyword == "dragmodel"){
        int x;
        list >> x;
        controller.set_dragmodel(x);
      }
      else if(keyword == "cd"){
        double x;
        list >> x;
        controller.set_Cd(x);
      }
      else if(keyword == "cdheight"){
        double x;
        list >> x;
        controller.set_cdheight(x);
      }
      else if(keyword == "visu_dt"){
        double dt;
        list >> dt;
        controller.set_visu_dt(dt);
      }
      else if(keyword == "visu_time"){
        double time;
        list >> time;
        controller.set_visu_time(time);
      }
      else if(keyword == "slice_dt"){
        double dt;
        list >> dt;
        controller.set_slice_dt(dt);
      }
      else if(keyword == "slice_time"){
        double time;
        list >> time;
        controller.set_slice_time(time);
      }
      else if(keyword == "uref"){
        double x;
        list >> x;
        controller.set_uref(x);
      }
      else if(keyword == "stat_time"){
        double x;
        list >> x;
        controller.set_stat_time(x);
      }
      else if(keyword == "vic"){
        double x1, x2, y1, y2, z1, z2, valx, valy, valz;
        list >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> valx >> valy >> valz;
        controller.set_vic(Vec3d(x1,y1,z1), Vec3d(x2,y2,z2), Vec3d(valx, valy, valz));  
      }
      else{
        if(mympi.master()){
          cerr << "Unknown input item ``" << keyword << "''\n";
        }
        exit(1);
      }
    } // end if keyword and control variable list
  } // end while
}

void Initialiser::init(MyMPI &mympi, Controller &controller, Grid &grid) const {

  // Sanity check
  bool error = false;
  for(int i=0; i<3; i++){
    if(grid_types[i] == 0) error = true;
  }
  if(error){
    if(mympi.rank() == 0) cerr << "Initialisation: Grid type error." << endl;
    exit(1);
  }
  for(int i=0; i<3; i++){
    if(mpi_blocks[i] <= 0) error = true;    
  }
  if(error){
    if(mympi.rank() == 0) cerr << "Initialisation: MPI blocks should be a positive integer." << endl;
    exit(1);
  }
  if(controller.nu() < 0){
    if(mympi.rank() == 0) cerr << "Initialisation: Did not specify a proper kinematic viscosity." << endl;
    exit(1);
  }
  if(controller.get_case() < 0){
    if(mympi.rank() == 0) cerr << "Initialisation: Error in specified case." << endl;
    exit(1);
  }
  if(controller.get_case() == CASE_URBAN){
    if(controller.utau() < 0){
      if(mympi.rank() == 0) cerr << "Initialisation: Did not specify friction velocity u_tau." << endl;
      exit(1);
    }
    if(controller.z0() < 0){
      if(mympi.rank() == 0) cerr << "Initialisation: Did not specify z0." << endl;
      exit(1);
    }
  }
  if(!controller.adaptive()){
    if(controller.dt() < 0){
      if(mympi.rank() == 0)
        cerr << "Initialisation: Did not specify a time step properly in a constant time stepping simulation." << endl;
      exit(1);
    }
  }
  if(controller.totaltime()<0){
    if(mympi.rank() == 0) cerr << "Initialisation: Did not specify a simulation total time properly." << endl;
    exit(1);
  }

  // Set controller
  // Pressure boundary conditions
  if((controller.get_case() == CASE_URBAN) || (controller.get_case() == CASE_TBL)){

    // Left pressure BC
    //controller.set_left_pbc(0, DIRICHLET);
    controller.set_left_pbc(0, NEUMANN  );
    controller.set_left_pbc(1, PERIODIC );
    controller.set_left_pbc(2, NEUMANN  );

    // Right pressure BC
    //controller.set_right_pbc(0, NEUMANN );
    controller.set_right_pbc(0, DIRICHLET);
    controller.set_right_pbc(1, PERIODIC );
    controller.set_right_pbc(2, NEUMANN  );

    // Left scalar BC
    controller.set_left_sbc(0, NEUMANN  );
    controller.set_left_sbc(1, PERIODIC );
    controller.set_left_sbc(2, NEUMANN  );

    // Right scalar BC
    controller.set_right_sbc(0, NEUMANN );
    controller.set_right_sbc(1, PERIODIC);
    controller.set_right_sbc(2, NEUMANN );
  }

  if(controller.get_case() == CASE_CYLINDER){

    // Left pressure BC
    controller.set_left_pbc(0, NEUMANN  );
    controller.set_left_pbc(1, PERIODIC );
    controller.set_left_pbc(2, NEUMANN  );

    // Right pressure BC
    controller.set_right_pbc(0, DIRICHLET);
    controller.set_right_pbc(1, PERIODIC );
    controller.set_right_pbc(2, NEUMANN  );

    // Left scalar BC
    controller.set_left_sbc(0, NEUMANN  );
    controller.set_left_sbc(1, PERIODIC );
    controller.set_left_sbc(2, NEUMANN  );

    // Right scalar BC
    controller.set_right_sbc(0, NEUMANN );
    controller.set_right_sbc(1, PERIODIC);
    controller.set_right_sbc(2, NEUMANN );
  }
  
  if(controller.get_case() == CASE_OPENCHANNEL){

    // Left pressure BC
    controller.set_left_pbc(0, PERIODIC );
    controller.set_left_pbc(1, PERIODIC );
    controller.set_left_pbc(2, NEUMANN  );

    // Right pressure BC
    controller.set_right_pbc(0, PERIODIC);
    controller.set_right_pbc(1, PERIODIC);
    controller.set_right_pbc(2, NEUMANN );

    // Left scalar BC
    controller.set_left_sbc(0, PERIODIC );
    controller.set_left_sbc(1, PERIODIC );
    controller.set_left_sbc(2, NEUMANN  );

    // Right scalar BC
    controller.set_right_sbc(0, PERIODIC);
    controller.set_right_sbc(1, PERIODIC);
    controller.set_right_sbc(2, NEUMANN );
  }

  // Make new directories 
  mkdir(  HDF_DIR, 0777);
  mkdir( GRID_DIR, 0777);
  mkdir( DATA_DIR, 0777);
  mkdir( SAVE_DIR, 0777);
  mkdir( STAT_DIR, 0777);
  mkdir(SLICE_DIR, 0777);
  mkdir(  VTK_DIR, 0777);

  // Initialise objects
  mympi.init(mpi_blocks);
  controller.init(mympi);   // init after BC are set
  grid.init(mympi);

  // Probes
  for(int i=0; i<probe_index.size(); i++){
    controller.set_probe(probe_index[i], probe_coords[i], grid, mympi);
  }

  // Slices
  for(int i=0; i<slice_index.size(); i++){
    controller.set_slice(slice_index[i], slice_direction[i], slice_coords[i], grid, mympi);
  }

}
