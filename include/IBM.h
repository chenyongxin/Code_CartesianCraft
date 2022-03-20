// IBM.h
// Immersed boundary method in 3D.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef IBM_H
#define IBM_H

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "Vec.h"
#include "Grid.h"
#include "Array.h"
#include "Field.h"
#include "MyMPI.h"
#include "parameter.h"
#include "Controller.h"
#include "VectorField.h"

using std::min;
using std::max;
using std::endl;
using std::cerr;
using std::string;
using std::ifstream;
using std::istringstream;

class IBM{

  private:

    // General variables.
    // Index [0,1,2] for vector and index 3 for scalar components
    Fieldd       _distance[4];     // signed distance function
    VectorFieldd _normal  [4];     // normalised normal vector
    Fieldb       _mask    [4];     // true for first grid point

    // Case specific variables
    // Cuboids
    int ncuboids;                  // number of cuboids
    Array2d xc;                    // centre coordinate, shape(xc) = [ncuboids, 3]
    Array2d dh;                    // cuboid side length, shape(dh) = [ncuboids, 3]

    // TODO: urban roughness

  public:

    IBM(){}
    ~IBM(){}

    // Initialise fields, read geom file and compute.
    void init(const Grid &, const Controller &, const MyMPI &);

    // Apply immersed boundary condition for velocity with no-slip BC.
    void ibc (VectorFieldd &v) const;               

    // Get signed distance.
    double distance(int d, int i, int j, int k) const {return _distance[d](i,j,k);  }

    // Get unit normal vector.
    Vec3d  normal  (int d, int i, int j, int k) const {return _normal[d].get(i,j,k);}

    // Get mask.
    bool   mask    (int d, int i, int j, int k) const {return _mask[d](i,j,k);      }

};

// Get distance and normal directions for an array of aligned cuboids.
void distance_normal_aligned_cuboids(int n, const Grid &grid, const Array2d &xc,
                                     const Array2d &dh, Fieldd *distance, VectorFieldd *normal);

// Get distance and normal directions for an aligned cuboid. 
void distance_normal_aligned_cuboid(const Vec3d &pt, const Vec3d &xc, const Vec3d &dh,
                                    double &distance, Vec3d &normal);

#endif // IBM_H
