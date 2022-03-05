// Psolver.h
// Poisson solver with HYPRE
// Code_CartesianCraft (c) Yongxin Chen

#ifndef PSOLVER_H
#define PSOLVER_H

#include <iostream>
#include <fstream>

#include "Vec.h"
#include "Grid.h"
#include "MyMPI.h"
#include "Field.h"
#include "parameter.h"
#include "Controller.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_sstruct_ls.h"

using std::ios;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;

class Psolver{

  private:

    // HYPRE objects
    HYPRE_SStructGrid     grid;
    HYPRE_SStructGraph    graph;
    HYPRE_SStructStencil  stencil;
    HYPRE_SStructMatrix   A;
    HYPRE_SStructVector   b;
    HYPRE_SStructVector   x;
    HYPRE_ParCSRMatrix    parA;
    HYPRE_ParVector       parb;
    HYPRE_ParVector       parx;
    HYPRE_Solver          solver;
    HYPRE_Solver          precond;
    HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};
    int                   object_type = HYPRE_PARCSR;

    // Only one part and one variable for HYPRE solver.
    int nparts    = 1;
    int nvars     = 1;
    int part      = 0;
    int var       = 0;

    // View
    int ilower[3];
    int iupper[3];

    // Data for passing rhs and collecting solution.
    double *data;

    // Flag, enquery if the object was initialised.
    bool initialised = false;

  public:
    Psolver();
    ~Psolver();

    // Initialise a Poisson solver.
    void init(const MyMPI &, const Grid &, const Controller &);

    // Solve linear equation Ap=rhs.
    void solve(const Fieldd &rhs, Fieldd &p, const MyMPI &); 
};

#endif // PSOLVER_H
