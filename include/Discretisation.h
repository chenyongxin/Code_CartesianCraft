// Discretisation.h
// Namespace to discretise the governing equations.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef DISCRETISATION_H
#define DISCRETISATION_H

#include <iostream>

#include "Vec.h"
#include "Grid.h"
#include "MyMPI.h"
#include "Field.h"
#include "parameter.h"
#include "Controller.h"
#include "VectorField.h"

namespace Discretisation{

  // Project pressure `p' to velocity field `v'.
  void project(VectorFieldd &v, const Fieldd &p, const Grid &grid, const Controller &controller);

  // Compute velocity divergence.
  void divergence(Fieldd &f, const VectorFieldd &v, const Grid &grid);

  // Compute velocity gradient tensor. Both staggered and central tensor are computed.
  void gradtensor(VectorFieldd *S_stag, VectorFieldd *S_cen, const VectorFieldd &v, const Grid &grid);

  // Compute convective-diffisive term and saved in VectorField a.
  void conv_diff(VectorFieldd &a, const VectorFieldd &v, const Grid &grid, const double nu, const Fieldd &ed, const VectorFieldd *uij);

  // Adam-Bashforth integration: v = v + (1.5r - 0.5h + bforce)*dt
  void AB(VectorFieldd &v, const VectorFieldd &r, const VectorFieldd &h, const Vec3d &bforce, const double dt);

};
#endif // DISCRETISATION