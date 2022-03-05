// Stat.cc
// Class holds time average statistics operations
// Code_CartesianCraft (c) Yongxin Chen

#include "Stat.h"
#ifdef USE_HDF
#include <mpi.h>
#include "hdf5.h"
#endif

void Stat::init(const Grid &grid){

  // Init fields
  u.allocate(grid.n());
  p.allocate(grid.n());
  for(int d=0; d<3; d++){
    uiui[d].allocate(grid.n());
    uiuj[d].allocate(grid.n());
  }
}

void Stat::do_statistics(const VectorFieldd &u, const Fieldd &p){

  for(int d=0; d<3; d++){
    mean1(this->u[d], u[d], samples);  
    mean2(uiui[d], u[d], u[d], samples);
  }
  mean2(uiuj[0], u[0], u[1], samples);
  mean2(uiuj[1], u[0], u[2], samples); 
  mean2(uiuj[2], u[1], u[2], samples);

  mean1(this->p, p, samples);

  samples++;
}

void Stat::write(const Grid &grid, const MyMPI &mympi) const {
   
  vector<Fieldd> fields;
  fields.reserve(10);         // reserve a big enough place first
  fields.push_back(uiui[0]);
  fields.push_back(uiui[1]);
  fields.push_back(uiui[2]);
  fields.push_back(uiuj[0]);
  fields.push_back(uiuj[1]);
  fields.push_back(uiuj[2]);
  
  vector<string> fieldnames;
  fieldnames.push_back(string("u1u1"));
  fieldnames.push_back(string("u2u2"));
  fieldnames.push_back(string("u3u3"));
  fieldnames.push_back(string("u1u2"));
  fieldnames.push_back(string("u1u3"));
  fieldnames.push_back(string("u2u3"));

  // Do substraction
  for(int i=0; i<uiui[0].data().size(); i++){
    fields[0].data().data()[i] -= u[0].data().data()[i]*u[0].data().data()[i];
    fields[1].data().data()[i] -= u[1].data().data()[i]*u[1].data().data()[i];
    fields[2].data().data()[i] -= u[2].data().data()[i]*u[2].data().data()[i];
    fields[3].data().data()[i] -= u[0].data().data()[i]*u[1].data().data()[i];
    fields[4].data().data()[i] -= u[0].data().data()[i]*u[2].data().data()[i];
    fields[5].data().data()[i] -= u[1].data().data()[i]*u[2].data().data()[i];
  }

  // Visualise it
#ifdef USE_HDF
  H5fluid_fields(u, p, fields, fieldnames, grid, mympi, (string(HDF_DIR)+string("/")).c_str(), "stat");
#else
  myvtk.fluid_fields(u, p, fields, fieldnames, grid, mympi, (string(STAT_DIR)+string("/")).c_str(), "stat");
#endif
}

void Stat::save_hdf(const MyMPI &mympi, const Grid &grid) const {
#ifdef USE_HDF  
  // Collectively write fields
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate( (string(HDF_DIR)+string("/")+string("Stat_save.h5")).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  H5Pclose(fapl);
  
  H5write_field(file_id, "u0"  , mympi, grid, u[0]   );
  H5write_field(file_id, "u1"  , mympi, grid, u[1]   );
  H5write_field(file_id, "u2"  , mympi, grid, u[2]   );

  H5write_field(file_id, "p"   , mympi, grid, p      );
  
  H5write_field(file_id, "u1u1", mympi, grid, uiui[0]);
  H5write_field(file_id, "u2u2", mympi, grid, uiui[1]);
  H5write_field(file_id, "u3u3", mympi, grid, uiui[2]);
  
  H5write_field(file_id, "u1u2", mympi, grid, uiuj[0]);
  H5write_field(file_id, "u1u3", mympi, grid, uiuj[1]);
  H5write_field(file_id, "u2u3", mympi, grid, uiuj[2]);

  H5Fclose(file_id);

  // Write step info, only once
  if(mympi.master()){
    ofstream file;
    file.open( (string(SAVE_DIR)+string("/")+string("Stat")).c_str(), ios::binary );
    file.write((char*)&samples, sizeof(int));
    file.close();
  }
#endif
}

void Stat::save(const MyMPI &mympi) const {
  ofstream file;
  char name[100];
  sprintf(name, "%s/Stat_x%dx%dx%d", SAVE_DIR, mympi.coords(0), mympi.coords(1), mympi.coords(2));  
  file.open(name, ios::binary);
  int steps = samples;
  file.write((char*)&steps, sizeof(int));
  u.save(file);
  p.save(file);
  for(int i = 0; i < 3; i++){
    uiui[i].save(file);
  }
  for(int i = 0; i < 3; i++){
    uiuj[i].save(file);
  }
  file.close();
}

void Stat::read_hdf(const MyMPI &mympi, const Grid &grid){
#ifdef USE_HDF
  // Collectively read fields
  hid_t file_id = H5Fopen( (string(HDF_DIR)+string("/")+string("Stat_save.h5")).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  H5read_field(file_id, "u0"  , mympi, grid, u[0]   );
  H5read_field(file_id, "u1"  , mympi, grid, u[1]   );
  H5read_field(file_id, "u2"  , mympi, grid, u[2]   );

  H5read_field(file_id, "p"   , mympi, grid, p      );

  H5read_field(file_id, "u1u1", mympi, grid, uiui[0]);
  H5read_field(file_id, "u2u2", mympi, grid, uiui[1]);
  H5read_field(file_id, "u3u3", mympi, grid, uiui[2]);

  H5read_field(file_id, "u1u2", mympi, grid, uiuj[0]);
  H5read_field(file_id, "u1u3", mympi, grid, uiuj[1]);
  H5read_field(file_id, "u2u3", mympi, grid, uiuj[2]);

  H5Fclose(file_id);

  // Read step info by all 
  ifstream file;
  file.open( (string(SAVE_DIR)+string("/")+string("Stat")).c_str(), ios::binary );
  int steps;
  file.read((char*)&steps, sizeof(int));
  samples = steps;
  file.close();
#endif
}

void Stat::read(const MyMPI &mympi){
  ifstream file;
  char name[100];
  sprintf(name, "%s/Stat_x%dx%dx%d", SAVE_DIR, mympi.coords(0), mympi.coords(1), mympi.coords(2));  
  file.open(name, ios::binary);
  int steps;
  file.read((char*)&steps, sizeof(int));
  samples = steps;
  for(int i = 0; i < 3; i++){
    uiui[i].read(file);
  } 
  for(int i = 0; i < 3; i++){
    uiuj[i].read(file);
  }
  file.close();
}

void mean1(Fieldd &mf, const Fieldd &f, int n){
  int size = mf.data().size();
  for(int i=0; i<size; i++){
    mf.data().data()[i] = (mf.data().data()[i]*n + f.data().data()[i])/(n+1);
  }
}

void mean2(Fieldd &mf, const Fieldd &f1, const Fieldd &f2, int n){
  int size = mf.data().size();
  for(int i=0; i<size; i++){
    mf.data().data()[i] = (mf.data().data()[i]*n + f1.data().data()[i] * f2.data().data()[i])/(n+1);
  }
}
