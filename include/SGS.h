// SGS.h
// Sub-grid scale model for LES.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef SGS_H
#define SGS_H

#include <cmath>

#include "Grid.h"
#include "Field.h"
#include "parameter.h"
#include "Controller.h"
#include "VectorField.h"

namespace SGS{

  // Compute subgrid-scale stress.
  // TODO: Add multiple SGS options.
  // ed: SGS eddy viscosity.
  // v: velocity
  // S: velocity gradient tensor.
  void compute_SGS(Fieldd &ed, const VectorFieldd &v, const VectorFieldd *S, const Grid &grid, const Controller &controller);

  // Vreman SGS model.
  // An eddy-viscosity subgrid-scale model for turbulent shear flow:
  // Algebraic theory and applications. PoF, 2004
  void Vreman(Fieldd &ed, const VectorFieldd *S, const Grid &grid);

}

#endif // SGS_H
