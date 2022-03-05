// H5IO.cc
// Write/Read field data to/from a HDF5 file collectively.
// Code_CartesianCraft (c) Yongxin Chen

#include "H5IO.h"

#ifdef USE_HDF
using namespace std;

// Write a double precision field
void H5write_field(hid_t file_id, const char *dataset_name,
                   const MyMPI &mympi, const Grid &grid, const Fieldd &f){

  hsize_t count[3] = {1,1,1};      // 1 piece to write
  hsize_t gn[3], n[3], offset[3];  // global, local dimensions and offset.
  int rbound[3] = {0,0,0};         // right boundary increment

  // Get data
  for(int d=0; d<3; d++){
    gn[d] = grid.gn(d)+1;          // vertex centre
    n[d]  = grid.n(d);             // number of cells
    if(mympi.rboundary(d)){
      n[d] += 1;
      rbound[d] = 1;
    }
    offset[d] = grid.piece_start(d); 
  }
  int size = n[0]*n[1]*n[2];

  // Preparation, transpose
  double *buffer = new double[size];
  int m = 0;
  for(int i=f.start(0); i<f.end(0)+rbound[0]; i++){
    for(int j=f.start(1); j<f.end(1)+rbound[1]; j++){
      for(int k=f.start(2); k<f.end(2)+rbound[2]; k++){
        buffer[m] = f(i, j, k);
        m++;
      } 
    }
  }

  const int NDIMS = 3; 
  hid_t fspace = H5Screate_simple(NDIMS, gn, NULL);
  hid_t dset = H5Dcreate(file_id, dataset_name, H5T_NATIVE_DOUBLE, fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Write data
  hid_t mspace = H5Screate_simple(NDIMS, n, NULL);
  H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, NULL, count, n);
  hid_t dcpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dcpl, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, dcpl, buffer);

  // Close things
  H5Sclose(mspace);
  H5Dclose(dset);
  H5Sclose(fspace);
  H5Pclose(dcpl);

  delete [] buffer;
}

// Read a double precision field 
void H5read_field(hid_t file_id, const char *dataset_name,
                  const MyMPI &mympi, const Grid &grid, Fieldd &f){

  hsize_t count[3] = {1,1,1};      // 1 piece to read
  hsize_t gn[3], n[3], offset[3];  // global, local dimensions and offset.
  int rbound[3] = {0,0,0};         // right boundary increment

  // Get data
  for(int d=0; d<3; d++){
    gn[d] = grid.gn(d)+1;          // vertex centre
    n[d]  = grid.n(d);             // number of cells
    if(mympi.rboundary(d)){
      n[d] += 1;
      rbound[d] = 1;
    }
    offset[d] = grid.piece_start(d); 
  }
  int size = n[0]*n[1]*n[2];
  
  // Retrieve data
  double *buffer = new double[size];
  hid_t dset = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
  hid_t fspace = H5Screate_simple(3, gn, NULL);
  hid_t mspace = H5Screate_simple(3, n, NULL);
  H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, NULL, count, n);
  H5Dread(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, buffer);

  // Transpose
  int m = 0;
  for(int i=f.start(0); i<f.end(0)+rbound[0]; i++){
    for(int j=f.start(1); j<f.end(1)+rbound[1]; j++){
      for(int k=f.start(2); k<f.end(2)+rbound[2]; k++){
        f(i, j, k) = buffer[m];
        m++;
      } 
    }
  }

  // Close things
  H5Sclose(mspace);
  H5Dclose(dset);
  H5Sclose(fspace);

  delete [] buffer;
}

// Write a single precision field
static void H5write_field_float(hid_t file_id, const char *dataset_name,
                   const MyMPI &mympi, const Grid &grid, const Fieldd &f){
  
  hsize_t count[3] = {1,1,1};
  hsize_t gn[3], n[3], offset[3];  // global, local dimensions and offset.
  int rbound[3] = {0,0,0};         // right boundary increment

  // Get data
  for(int d=0; d<3; d++){
    gn[d] = grid.gn(d)+1;          // vertex centre
    n[d]  = grid.n(d);             // number of cells
    if(mympi.rboundary(d)){
      n[d] += 1;
      rbound[d] = 1;
    }
    offset[d] = grid.piece_start(d); 
  }
  int size = n[0]*n[1]*n[2];

  // Preparation
  float *buffer = new float[size];
  int m = 0;
  for(int i=f.start(0); i<f.end(0)+rbound[0]; i++){
    for(int j=f.start(1); j<f.end(1)+rbound[1]; j++){
      for(int k=f.start(2); k<f.end(2)+rbound[2]; k++){
        buffer[m] = f(i, j, k);
        m++;
      } 
    }
  }

  const int NDIMS = 3; 
  hid_t fspace = H5Screate_simple(NDIMS, gn, NULL);
  hid_t dset = H5Dcreate(file_id, dataset_name, H5T_NATIVE_FLOAT, fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Write data
  hid_t mspace = H5Screate_simple(NDIMS, n, NULL);
  H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, NULL, count, n);

  hid_t dcpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dcpl, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset, H5T_NATIVE_FLOAT, mspace, fspace, dcpl, buffer);

  // Close things
  H5Sclose(mspace);
  H5Dclose(dset);
  H5Sclose(fspace);
  H5Pclose(dcpl);

  delete [] buffer;
}



static void H5write_grid(const char *filename, const Grid &grid, const MyMPI &mympi){

  // 1. Open file collectively
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  H5Pclose(fapl);

  // 2. Get data and write
  Fieldd xyz(grid.n());
  const char *coords[3] = {"x", "y", "z"}; 
  for(int d=0; d<3; d++){
    for(int k=xyz.ks(); k<xyz.ke()+1; k++){
      for(int j=xyz.js(); j<xyz.je()+1; j++){
        for(int i=xyz.is(); i<xyz.ie()+1; i++){
          int ijk[3] = {i, j, k};
          xyz(i,j,k) = grid.x(d, ijk[d]);
        }
      }
    }
    H5write_field_float(file_id, coords[d], mympi, grid, xyz);
  }

  // 3. Close file
  H5Fclose(file_id);
}


// A wrapper of writing single precision fluid fields
void H5fluid_fields(const VectorFieldd   &u                 ,
                    const Fieldd         &p                 ,
                    const vector<Fieldd> &other             ,
                    const vector<string> &fieldnames        ,
                    const Grid           &grid              ,
                    const MyMPI          &mympi             ,
                    const char           *prefix            ,
                    const char           *filename          ){

  // Retrieve number of fields
  int nfields = other.size();

  // Make file name 
  char name[200];
  sprintf(name, "%s%s.h5", prefix, filename);

  // Part 1: Output data
  // 1. Open file collectively
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  H5Pclose(fapl);

  // 2. Write data collectively
  H5write_field_float(file_id, "u", mympi, grid, u[0]);
  H5write_field_float(file_id, "v", mympi, grid, u[1]);
  H5write_field_float(file_id, "w", mympi, grid, u[2]);
  H5write_field_float(file_id, "p", mympi, grid, p   );
  for(int i=0; i<nfields; i++){ 
    H5write_field_float(file_id, fieldnames[i].c_str(), mympi, grid, other[i]);
  }

  // 3. Close file
  H5Fclose(file_id);

  // Part 2: Output grid only once
  static bool made_grid = false;
  if(!made_grid){
    sprintf(name, "%sgrid.h5", prefix);
    H5write_grid(name, grid, mympi);
    made_grid = true;
  }

  // Part 3: Write xdmf file
  if(mympi.master()){
    char buffer[200];
    sprintf(name, "%sf%s.xdmf", prefix, filename);           // {prefix}/f{filename}.{ext} with a `f` in the middle
    ofstream xdmf;
    xdmf.open(name);
    xdmf << "<?xml version=\"1.0\" ?>" << endl; 
    xdmf << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << endl;
    xdmf << "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">" << endl;
    xdmf << "<Domain>" << endl;
    xdmf << "<Grid Name=\"Structured Grid\" GridType=\"Uniform\">" << endl; 
    sprintf(buffer, "<Topology TopologyType=\"3DSMesh\" Dimensions=\"%d %d %d\"/>", grid.gn(0)+1, grid.gn(1)+1, grid.gn(2)+1);
    xdmf << buffer << endl;
    xdmf << "<Geometry name=\"geo\" Type=\"X_Y_Z\">" << endl;
    char coords[3] = {'x', 'y', 'z'};
    for(int d=0; d<3; d++){
      sprintf(buffer, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\ngrid.h5:/%c\n</DataItem>",
          grid.gn(0)+1, grid.gn(1)+1, grid.gn(2)+1, coords[d]);
      xdmf << buffer << endl;
    }
    xdmf << "</Geometry>" << endl;
    
    // Velocity vector
    xdmf << "<Attribute Name=\"u\" AttributeType=\"Vector\" Center=\"Node\">" << endl;
    sprintf(buffer, "<DataItem ItemType=\"Function\" Dimensions=\"%d %d %d 3\" Function=\"JOIN($0 , $1 , $2)\">", 
        grid.gn(0)+1, grid.gn(1)+1, grid.gn(2)+1);
    xdmf << buffer << endl;
    
    // u
    sprintf(buffer, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n%s.h5:/%s\n</DataItem>",
        grid.gn(0)+1, grid.gn(1)+1, grid.gn(2)+1, filename, "u");
    xdmf << buffer << endl;
    
    // v
    sprintf(buffer, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n%s.h5:/%s\n</DataItem>",
        grid.gn(0)+1, grid.gn(1)+1, grid.gn(2)+1, filename, "v");
    xdmf << buffer << endl;
    
    // w
    sprintf(buffer, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n%s.h5:/%s\n</DataItem>",
        grid.gn(0)+1, grid.gn(1)+1, grid.gn(2)+1, filename, "w");
    xdmf << buffer << endl;
    
    xdmf << "</DataItem>" << endl << "</Attribute>" << endl;
    
    // Pressure
    xdmf << "<Attribute Name=\"p\" Center=\"Node\">" << endl;
    sprintf(buffer, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n%s.h5:/%s\n</DataItem>",
        grid.gn(0)+1, grid.gn(1)+1, grid.gn(2)+1, filename, "p");
    xdmf << buffer << endl;
    xdmf << "</Attribute>" << endl;

    // Other fields
    for(int i=0; i<nfields; i++){ 
      xdmf << "<Attribute Name=\""<< fieldnames[i] << "\" Center=\"Node\">" << endl;
      sprintf(buffer, "<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n%s.h5:/%s\n</DataItem>",
          grid.gn(0)+1, grid.gn(1)+1, grid.gn(2)+1, filename, fieldnames[i].c_str());
      xdmf << buffer << endl;
      xdmf << "</Attribute>" << endl;
    }

    // Ending
    xdmf << "</Grid>" << endl;
    xdmf << "</Domain>" << endl;
    xdmf << "</Xdmf>" << endl;

    // Close file
    xdmf.close();
  }

}

#endif // USE_HDF
