// BC.h
// Boundary condition class.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef BC_H
#define BC_H

#include <cmath>
#include <random>
#include <string>
#include <sstream>
#include <fstream>

#include "Vec.h"
#include "Grid.h"
#include "Array.h"
#include "Field.h"
#include "MyMPI.h"
#include "parameter.h"
#include "Controller.h"
#include "VectorField.h"

using std::ios;
using std::string;
using std::ifstream;
using std::istringstream;

class BC{

  private:

    // View for MPI data exchange
    Vec3i left_view_set_start    [3];      // view for logical left set ghost cell start and end
    Vec3i left_view_set_end      [3];      // exterior

    Vec3i left_view_get_start    [3];      // view for logical left get ghost cell start and end
    Vec3i left_view_get_end      [3];      // interior

    Vec3i right_view_set_start   [3];      // view for logical right set ghost cell start and end
    Vec3i right_view_set_end     [3];      // exterior

    Vec3i right_view_get_start   [3];      // view for logical right set ghost cell start and end
    Vec3i right_view_get_end     [3];      // interior

    // View for slices: -1 offset ghost cell layer
    Vec3i  left_view_slice_start [3];
    Vec3i  left_view_slice_end   [3];
    Vec3i right_view_slice_start [3];
    Vec3i right_view_slice_end   [3];

    // MPI rank for different BC types
    Vec3i rank_left_vbc;
    Vec3i rank_right_vbc;
    Vec3i rank_left_pbc;
    Vec3i rank_right_pbc;
    Vec3i rank_left_sbc;
    Vec3i rank_right_sbc;

    // Slice values for applying boundary conditions: xxx_slice[direction][component]
    Array3d  left_slice[3][3];
    Array3d right_slice[3][3];

    // Projected area. A slice that includes ghost cells.
    Array3d area[3];

    // Turbulent inflow generation
    // Profiles along the height
    Array1d  umean,  vmean,  wmean;       // local mean velocity profile, *mean[nz]
    Array1d gumean, gvmean, gwmean;       // global mean velocity profile
    Array1d  uu, vv, ww;                  // global normal Reynolds stress, xx[nz]
    Array1d  uv, uw, vw;                  // global shear Reynolds stress, xy[nz]
    Array1d  Lt, Lx, Ly, Lz;              // global integral length scales, Lx[nz] and Langragian time scale
    Array1d a11,a21,a22;                  // global Lund matrix
    Array1d a31,a32,a33;                  // global Lund matrix
    Array1i  ny, nz;                      // global integral length scale in grid cells
    // 2D slice
    Array2d random_u, random_v, random_w; // global random slice for u,v,w
    Array2d bfilty, bfiltz;               // global filter in y and z, bij = bi*bj
    Array2d psiu, psiv, psiw;             // global Psi slice
    Array2d psiu_old, psiv_old, psiw_old; // global Psi slice in the last step
    // Fluctuation slice
    Array2d  uf,  vf,  wf;                // local fluctuation for u,v,w, same size with inner inlet plane
    Array2d guf, gvf, gwf;                // global velocity fluctuation planes

    // Random number generator
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> distribution{0,1}; // 0 mean, 1 variance distribution

    // Set view for MPI data exchange regions and boundary slices.
    void set_view(const Grid &);

    // Read and assign profile along the whole height for turbulent inflow generator.
    void read_turbulent_inflow_file(Array1d &var, const Grid &grid, const char *filename);

    // Initialise turbulent inflow generation related variables.
    void init_turbulent_inflow_generation(const Grid &);

    // Generate random number to random slices and update psi
    void random_gen();

  public:

    BC(){};
    ~BC(){}

    // Initialise a BC object and setup ranks for specific BC types.
    void init     (const MyMPI &, const Controller &, const Grid &);

    // Apply velocity bounday condition.
    void apply_vbc(VectorFieldd &v, const MyMPI &, const Controller &, const Grid &);

    // Apply pressure boundary condition.
    void apply_pbc(      Fieldd &f, const MyMPI &, const Controller &) const;

    // Apply boundary condition for scalars.
    void apply_sbc(      Fieldd &f, const MyMPI &, const Controller &) const;

    // Apply velocity initial condition by using boundary condition.
    void apply_vic(VectorFieldd &v, const Controller &, const Grid &);

    // Generate new velocity fluctuation.
    void fluctuation_gen(const MyMPI &, const Controller &, const Grid &);

    // Halo cells swap with specific left and right ranks
    template<class T> void swap(Field<T> &f, const MyMPI &, const Vec3i &left, const Vec3i &right) const;

};

// Function template for swapping ghost cell data.
template<class T>
void BC::swap(Field<T> &f, const MyMPI &mympi, const Vec3i &left, const Vec3i &right) const {

  // Get left halo; send to the left and receive from the right
  for(int d=0; d<3; d++){
    int size = (left_view_get_end[d] - left_view_get_start[d]).product();
    T *sendbuff = new T[size];
    T *recvbuff = new T[size];
    f.pack(left_view_get_start[d], left_view_get_end[d], sendbuff);
    mympi.sendrecv(sendbuff, recvbuff, size, mympi.left(d), mympi.right(d));
    if(right(d) >= 0){ // get data but only fill into the field as it has the right
      f.fill(right_view_set_start[d], right_view_set_end[d], recvbuff);
    }
    delete [] sendbuff;
    delete [] recvbuff;
  }

  // Get right halo; send to the right and recive from the left
  for(int d=0; d<3; d++){
    int size = (right_view_get_end[d] - right_view_get_start[d]).product();
    T *sendbuff = new T[size];
    T *recvbuff = new T[size];
    f.pack(right_view_get_start[d], right_view_get_end[d], sendbuff);
    mympi.sendrecv(sendbuff, recvbuff, size, mympi.right(d), mympi.left(d));
    if(left(d) >= 0){ // get data but only fill into the field as it has the left
      f.fill(left_view_set_start[d], left_view_set_end[d], recvbuff);
    }
    delete [] sendbuff;
    delete [] recvbuff;
  }
}

// Neumann BC
void Neumann_bc(Fieldd &f,
                const Vec3i &set_start, const Vec3i &set_end,
                const Vec3i &get_start, const Vec3i &get_end);

// Dirichlet BC for surface in between cell centres with value 0
void Dirichlet_bc(Fieldd &f,
                  const Vec3i &set_start, const Vec3i &set_end,
                  const Vec3i &get_start, const Vec3i &get_end);

// Dirichlet BC for a surface in between cell centres with a specified value
void Dirichlet_bc(Fieldd &f,
                  const Vec3i &set_start, const Vec3i &set_end,
                  const Vec3i &get_start, const Vec3i &get_end,
                  const double *data);

// Dirichlet BC for a surface with a specified value
void Dirichlet_bc(Fieldd &f,
                  const Vec3i &set_start, const Vec3i &set_end, const double *data);

// Dirichlet BC for a specified value on the slice
void Dirichlet_bc(Fieldd &f, const Vec3i &set_start, const Vec3i &set_end,
                  const double value);

#endif // BC_H
