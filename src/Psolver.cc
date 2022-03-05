// Psolver.cc
// Poisson solver with HYPRE
// Code_CartesianCraft (c) Yongxin Chen

#include "Psolver.h"

Psolver::Psolver(){}

Psolver::~Psolver(){
  if(initialised){
    HYPRE_ParCSRPCGDestroy(solver);
    HYPRE_BoomerAMGDestroy(precond);
  }
}

void Psolver::init(const MyMPI &mympi, const Grid &mygrid, const Controller &controller){

  // Dimension and stencils
  const int ndims     = 3;
  const int nstencils = 7;
  int stencil_indices[nstencils] = {0,1,2,3,4,5,6};
  int offsets[nstencils][ndims]  = {{ 0, 0, 0},
                                    {-1, 0, 0},
                                    { 0,-1, 0},
                                    { 0, 0,-1},
                                    { 1, 0, 0},
                                    { 0, 1, 0},
                                    { 0, 0, 1}};

  // Periodicity
  int periodic[3] = {0, 0, 0};

  // 0. Get control variables.
  // 0.1 Set view
  for(int d=0; d<3; d++){
    ilower[d] = mygrid.piece_start(d);
    iupper[d] = mygrid.piece_end(d)-1;
  }
  // 0.2 Set up stencils
  Fieldd left[3], right[3], diagonal;  // left, right off diagonal and diagonal matrix.
  diagonal.allocate(mygrid.n());
  for(int d=0; d<3; d++){
    left[d].allocate(mygrid.n());
    right[d].allocate(mygrid.n());
  }
  for(int d=0; d<3; d++){
    Vec3i o(0,0,0); o(d) = 1;    // offset
    for(int k=mygrid.start(2); k<mygrid.end(2); k++){
      for(int j=mygrid.start(1); j<mygrid.end(1); j++){
        for(int i=mygrid.start(0); i<mygrid.end(0); i++){
          Vec3i ijk(i,j,k);      // ijk
          Vec3i oijk = ijk+o;    // offseted ijk
          left[d](i,j,k)   = -mygrid.dxi(d, ijk(d))*mygrid.dxi2(d,  ijk(d))*mygrid.dv(i,j,k);
          right[d](i,j,k)  = -mygrid.dxi(d, ijk(d))*mygrid.dxi2(d, oijk(d))*mygrid.dv(i,j,k);
          diagonal(i,j,k) -= left[d](i,j,k) + right[d](i,j,k);
        }
      }
    }
  }
  // 0.3 Modify stencils for boundary conditions  
  // Apply boundary conditions to modify stencil
  for(int d=0; d<3; d++){

    // Periodicity
    // Periodicity in this direction?
    bool periodicity = (controller.left_pbc(d) == PERIODIC) || (controller.right_pbc(d)); // single partition check
    mympi.lor(periodicity);                                                               // over all the partitions
    if(periodicity) periodic[d] = mygrid.gn(d);
    // Left boundaries
    {
      Vec3i start = mygrid.start();
      Vec3i end   = mygrid.end()  ;
      end(d)      = start(d)+1;         // flatten to starting slice

      for(int k=start(2); k<end(2); k++){
        for(int j=start(1); j<end(1); j++){
          for(int i=start(0); i<end(0); i++){
            if(controller.left_pbc(d) == NEUMANN){
              diagonal(i,j,k) += left[d](i,j,k);
              left[d](i,j,k)   = 0;
            }
            else if(controller.left_pbc(d) == DIRICHLET){
              diagonal(i,j,k) -= left[d](i,j,k);
              left[d](i,j,k)   = 0;
            }
          }
        }
      }
    }
    // Right boundaries
    {
      Vec3i start = mygrid.start();
      Vec3i end   = mygrid.end()  ;
      start(d)    = end(d)-1;          // flatten to ending slice 

      for(int k=start(2); k<end(2); k++){
        for(int j=start(1); j<end(1); j++){
          for(int i=start(0); i<end(0); i++){
            if(controller.right_pbc(d) == NEUMANN){
              if(mympi.rboundary(d)){
                diagonal(i,j,k) += right[d](i,j,k);
                right[d](i,j,k)  = 0;
              }
            }
            else if(controller.right_pbc(d) == DIRICHLET){
              if(mympi.rboundary(d)){
                diagonal(i,j,k) -= right[d](i,j,k);
                right[d](i,j,k)  = 0;
              }
            }
          }
        }
      }
    }
  }

  // Flatten fields and take data.
  double *left1, *left2, *left3, *right1, *right2, *right3, *diag;
  int size = mygrid.size();
  left1  = new double[size];
  left2  = new double[size];
  left3  = new double[size];
  right1 = new double[size];
  right2 = new double[size];
  right3 = new double[size];
  diag   = new double[size];

  Vec3i start = mygrid.start();
  Vec3i end   = mygrid.end()  ;
  diagonal.pack(start, end, diag  );
  left[0].pack( start, end, left1 );
  left[1].pack( start, end, left2 );
  left[2].pack( start, end, left3 );
  right[0].pack(start, end, right1);
  right[1].pack(start, end, right2);
  right[2].pack(start, end, right3);

  // 0.4 Allocate data
  data = new double[size];
  for(int i=0; i<size; i++) data[i] = 0;

  // 1. Set up the 3D grid.
  HYPRE_SStructGridCreate(mympi.comm(), ndims, nparts, &grid);
  HYPRE_SStructGridSetExtents(grid, part, ilower, iupper);
  HYPRE_SStructGridSetVariables(grid, part, nvars, vartypes);
  HYPRE_SStructGridSetPeriodic(grid, part, periodic);
  HYPRE_SStructGridAssemble(grid);

  // 2. Define the discretisation stencils.
  HYPRE_SStructStencilCreate(ndims, nstencils, &stencil);
  for(int entry=0; entry<nstencils; entry++){
    HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], var);
  }

  // 3. Set up the Graph.
  HYPRE_SStructGraphCreate(mympi.comm(), grid, &graph);
  HYPRE_SStructGraphSetObjectType(graph, object_type);
  HYPRE_SStructGraphSetStencil(graph, part, var, stencil);
  HYPRE_SStructGraphAssemble(graph);

  // 4. Set up a SStruct Matrix.
  int    nvalues = size * nstencils;
  double *values = new double[nvalues];
  int m = 0;
  for(int i=0; i<nvalues; i+=nstencils){
    values[i  ] =   diag[m];
    values[i+1] =  left1[m];
    values[i+2] =  left2[m];
    values[i+3] =  left3[m];
    values[i+4] = right1[m];
    values[i+5] = right2[m];
    values[i+6] = right3[m];
    m++;
  }
  HYPRE_SStructMatrixCreate(mympi.comm(), graph, &A);
  HYPRE_SStructMatrixSetObjectType(A, object_type);
  HYPRE_SStructMatrixInitialize(A);
  HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, var, nstencils, 
      stencil_indices, values);
  HYPRE_SStructMatrixAssemble(A);

  // 5. Set up SStruct Vectors for b and x.
  HYPRE_SStructVectorCreate(mympi.comm(), grid, &b);
  HYPRE_SStructVectorCreate(mympi.comm(), grid, &x);
  HYPRE_SStructVectorSetObjectType(b, object_type);
  HYPRE_SStructVectorSetObjectType(x, object_type);
  HYPRE_SStructVectorInitialize(b);
  HYPRE_SStructVectorInitialize(x);
  HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, var, data); // initial value for b
  HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, var, data); // initial value for x
  HYPRE_SStructVectorAssemble(b);
  HYPRE_SStructVectorAssemble(x);

  // 6. Set up and use a solver.
  HYPRE_SStructMatrixGetObject(A, (void **) &parA);
  HYPRE_SStructVectorGetObject(b, (void **) &parb);
  HYPRE_SStructVectorGetObject(x, (void **) &parx);
  HYPRE_ParCSRPCGCreate(mympi.comm(), &solver);
  HYPRE_ParCSRPCGSetTol(solver, HYPRE_TOLERANCE);
  HYPRE_ParCSRPCGSetPrintLevel(solver, HYPRE_PCG_PRINT_LV);
  HYPRE_ParCSRPCGSetMaxIter(solver, HYPRE_MAX_ITERATION);
  HYPRE_BoomerAMGCreate(&precond);       // create the BoomerAMG solver for use as a preconditioner
  HYPRE_BoomerAMGSetMaxIter(precond, 1); // Set BoomerAMG parameters 
  HYPRE_BoomerAMGSetTol(precond, 0.0);
  HYPRE_BoomerAMGSetPrintLevel(precond, HYPRE_AMG_PRINT_LV); // print amg solution info 
  HYPRE_BoomerAMGSetRelaxType(precond, 6);                   // Sym G.S./Jacobi hybrid 
  HYPRE_BoomerAMGSetNumSweeps(precond, 1);
  HYPRE_ParCSRPCGSetPrecond(solver, HYPRE_BoomerAMGSolve,
      HYPRE_BoomerAMGSetup, precond);  // Set preconditioner
  HYPRE_ParCSRPCGSetup(solver, parA, parb, parx);

  // 7. Clean up memory
  delete [] diag;
  delete [] left1;
  delete [] left2;
  delete [] left3;
  delete [] right1;
  delete [] right2;
  delete [] right3;

  // 8. Update status
  initialised = true;
}


void Psolver::solve(const Fieldd &rhs, Fieldd &p, const MyMPI &mympi){

  // Iteration and residual
  int    num_iterations;
  double residual;

  // Solve
  rhs.pack(rhs.start(), rhs.end(), data);                              // get rhs data
  HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, var, data); // set rhs data to b
  HYPRE_SStructVectorAssemble(b);                                      // get the rhs ready  
  HYPRE_ParCSRPCGSolve(solver, parA, parb, parx);                      // solve the equations
  HYPRE_SStructVectorGather(x);                                        // gather data
  HYPRE_SStructVectorGetBoxValues(x, part, ilower, iupper, var, data); // get data from solution
  p.fill(p.start(), p.end(), data);                                    // fill data to field
  HYPRE_ParCSRPCGGetNumIterations(solver, &num_iterations);            // get number of iterations
  HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(solver, &residual);      // get residual

  // Output calculation result
  if(mympi.rank() == 0){
    ofstream file;
    file.open(MONITOR_HYPRE, ios::app);
    file << num_iterations << " " << " " << residual << endl;
    file.close();
  }

  // Check if diverges
  if((num_iterations >= HYPRE_MAX_ITERATION) && (residual >= HYPRE_TOLERANCE)){
    if(mympi.rank() == 0){
      cerr << "Poisson solver divergers!" << endl;
    }
    exit(1);
  }
}
