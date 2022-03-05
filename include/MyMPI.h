// MyMPI.h
// Customised MPI class.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef MYMPI_H
#define MYMPI_H

#include <mpi.h>
#include <iostream>
#include <typeinfo>

using std::cout;
using std::endl;
using std::cerr;

class MyMPI{

  private:
    int      _rank         = 0;                     // MPI rank
    int      _size         = 0;                     // MPI communicator world size
    MPI_Comm _comm;                                 // MPI_COMM_WORLD 
    int      _coords[3]    = { 0, 0, 0};            // MPI topo coordinates
    int      _left  [3]    = {-1,-1,-1};            // left neighbour rank
    int      _right [3]    = {-1,-1,-1};            // right neighbour rank
    int      _blocks[3]    = { 1, 1, 1};            // number of blocks in 3D

  public:

    // Constructor and deconstructor
    MyMPI(int &argc, char **argv);                  // MPI initialisation
    ~MyMPI();                                       // MPI finalisation

    // Initialise MPI with Cartesian topology.
    // Input is a dimension with 3 entries.
    void init(const int *dims);                     // called after MyMPI initialised

    // Enquiry
    int      rank     ()             const;         // get MPI rank
    int      size     ()             const;         // get number of processors
    MPI_Comm comm     ()             const;         // get MPI_COMM_WORLD 
    int      coords   (const int d)  const;         // get MPI coordinate in direction d
    int      left     (const int d)  const;         // get left neighbour's id in direction d
    int      right    (const int d)  const;         // get right neighbour's id in direction d
    int      blocks   (const int d)  const;         // get number of blocks in direction d
    bool     lboundary(const int d)  const;         // if left wall is the boundary
    bool     rboundary(const int d)  const;         // if right wall is the boundary
    bool     master   ()             const;         // if the current processor is the master processor

    // Collective operator
    template<class T> void max (T &) const;         // get the maximum of a specified value 
    template<class T> void min (T &) const;         // get the minimum of a specified value
    template<class T> void land(T &) const;         // logical and
    template<class T> void lor (T &) const;         // logical or
    template<class T> void sum (T &) const;         // sum 

    // Broadcast data from root to other processes
    template<class T> void bcast(T *data, int count, int root=0) const;

    // Send and receive data for halo cells
    template<class T> 
      void sendrecv(T *sendbuff, T *recvbuff, int count, int dest, int src) const; 

    // Wait until all processes reach the point
    void wait(){ MPI_Barrier(_comm); }
};


inline int      MyMPI::rank      ()            const {return _rank;                   }
inline int      MyMPI::size      ()            const {return _size;                   }
inline MPI_Comm MyMPI::comm      ()            const {return _comm;                   }
inline int      MyMPI::coords    (const int d) const {return _coords[d];              }
inline int      MyMPI::left      (const int d) const {return _left[d];                }
inline int      MyMPI::right     (const int d) const {return _right[d];               }
inline int      MyMPI::blocks    (const int d) const {return _blocks[d];              }
inline bool     MyMPI::lboundary (const int d) const {return _coords[d]==0;           }
inline bool     MyMPI::rboundary (const int d) const {return _coords[d]==_blocks[d]-1;}
inline bool     MyMPI::master    ()            const {return _rank==0;                }

// Function template for choosing a data type used in MPI
template<class T> MPI_Datatype choose_dtype(){
  if(typeid(T) == typeid(int)){
    return MPI_INT;
  }
  else if(typeid(T) == typeid(float)){
    return MPI_FLOAT;
  }
  else if(typeid(T) == typeid(double)){
    return MPI_DOUBLE;
  }
  else if(typeid(T) == typeid(bool)){
    return MPI_C_BOOL;
  }
  else{
    cerr << "Unsupported data type in MyMPI." << endl;
    exit(1);
  }
}

template<class T>
void MyMPI::bcast(T *data, int count, int root) const {
  MPI_Datatype dtype = choose_dtype<T>();
  MPI_Bcast(data, count, dtype, root, _comm);
}

template<class T>
void MyMPI::sendrecv(T *sendbuff, T *recvbuff, int count, int dest, int src) const {
  MPI_Request  requests[2];
  MPI_Status   statuses[2];
  MPI_Datatype dtype = choose_dtype<T>();
  MPI_Isend(sendbuff, count, dtype, dest, 1, _comm, &requests[0]);
  MPI_Irecv(recvbuff, count, dtype, src,  1, _comm, &requests[1]);
  MPI_Waitall(2, requests, statuses);
}

template<class T>
void MyMPI::max(T &b) const {
  T all, each = b;
  MPI_Datatype dtype = choose_dtype<T>();
  MPI_Allreduce(&each, &all, 1, dtype, MPI_MAX, _comm);
  b = all;
}

template<class T>
void MyMPI::min(T &b) const {
  T all, each = b;
  MPI_Datatype dtype = choose_dtype<T>();
  MPI_Allreduce(&each, &all, 1, dtype, MPI_MIN, _comm);
  b = all;
}

template<class T>
void MyMPI::land(T &b) const {
  T all, each = b;
  MPI_Datatype dtype = choose_dtype<T>();
  MPI_Allreduce(&each, &all, 1, dtype, MPI_LAND, _comm);
  b = all;
}

template<class T>
void MyMPI::lor(T &b) const {
  T all, each = b;
  MPI_Datatype dtype = choose_dtype<T>();
  MPI_Allreduce(&each, &all, 1, dtype, MPI_LOR, _comm);
  b = all;
}

template<class T>
void MyMPI::sum(T &b) const {
  T all, each = b;
  MPI_Datatype dtype = choose_dtype<T>();
  MPI_Allreduce(&each, &all, 1, dtype, MPI_SUM, _comm);
  b = all;
}

#endif // MYMPI_H
