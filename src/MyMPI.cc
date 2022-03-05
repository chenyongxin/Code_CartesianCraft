// MyMPI.cc
// Customised MPI class.
// Code_CartesianCraft (c) Yongxin Chen

#include "MyMPI.h"

MyMPI::MyMPI(int &argc, char **argv): _comm(MPI_COMM_WORLD) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(_comm, &_rank);
  MPI_Comm_size(_comm, &_size);
}

MyMPI::~MyMPI(){
  MPI_Finalize();
}

void MyMPI::init(const int *dims){
  int size = dims[0]*dims[1]*dims[2];
  if(size != _size){
    if(_rank == 0){
      cerr << "MPI processors not matched. Set " << size << " while " << _size << " given." << endl;
    }
    exit(1);
  }
  int periodic[3] = {1,1,1};
  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periodic, 1, &_comm);
  MPI_Cart_coords(_comm, _rank, 3, _coords);
  int l, r;
  for(int d=0; d<3; d++){
    _blocks[d] = dims[d];
    MPI_Cart_shift(_comm, d, 1, &l, &r);
    _left[d]  = l;
    _right[d] = r;
  }
}
