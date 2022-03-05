// BC.cc
// Boundary condition class.
// Code_CartesianCraft (c) Yongxin Chen

#include "BC.h"

void BC::set_view(const Grid &grid){

  // Get piece's number of cells
  Vec3i n = grid.n();

  // Left view set start
  for(int d=0; d<3; d++){
    Vec3i offset(0);
    left_view_set_start[d] = offset;
  }

  // Left view set end
  for(int d=0; d<3; d++){
    Vec3i offset(n+2*GC);
    offset(d) = GC;
    left_view_set_end[d] = offset;
  }

  // Left view get start
  for(int d=0; d<3; d++){
    Vec3i offset(0);
    offset(d) = GC;
    left_view_get_start[d] = offset;
  }

  // Left view get end
  for(int d=0; d<3; d++){
    Vec3i offset(n+2*GC);
    offset(d) = 2*GC;
    left_view_get_end[d] = offset;
  }

  // Right view set start
  for(int d=0; d<3; d++){
    Vec3i offset(0);
    offset(d) = n(d)+GC;
    right_view_set_start[d] = offset;
  }

  // Right view set end
  for(int d=0; d<3; d++){
    Vec3i offset(n+2*GC);
    right_view_set_end[d] = offset;
  }

  // Right view get start
  for(int d=0; d<3; d++){
    Vec3i offset(0);
    offset(d) = n(d);
    right_view_get_start[d] = offset;
  }

  // Right view get end
  for(int d=0; d<3; d++){
    Vec3i offset(n+2*GC);
    offset(d) = n(d)+GC;
    right_view_get_end[d] = offset;
  }

  // Left view slice start
  for(int d=0; d<3; d++){
    Vec3i offset(0);
    offset(d) = GC-1;
    left_view_slice_start[d] = offset;
  }

  // Left view slice end
  for(int d=0; d<3; d++){
    Vec3i offset(n+2*GC);
    offset(d) = GC;
    left_view_slice_end[d] = offset;
  }

  // Right view slice start
  for(int d=0; d<3; d++){
    Vec3i offset(0);
    offset(d) = n(d)+GC;
    right_view_slice_start[d] = offset;
  }

  // Right view slice end
  for(int d=0; d<3; d++){
    Vec3i offset(n+2*GC);
    offset(d) = n(d)+GC+1;
    right_view_slice_end[d] = offset;
  }

}

void BC::init(const MyMPI &mympi, const Controller &controller, const Grid &grid){

  // Set view for MPI data exchange and slices
  set_view(grid);

  // Init boundary slices and projected area slice
  Vec3i n = grid.n();
  for(int d=0; d<3; d++){
    Vec3i dims = n+2*GC;        // include ghost cells
    dims(d) = 1;                // flatten field to a slice in normal direction

    // Assistant velocity slice
    for(int e=0; e<3; e++){
      left_slice [d][e].allocate(dims(0), dims(1), dims(2));   // it is array not field
      right_slice[d][e].allocate(dims(0), dims(1), dims(2));   // it is the whole plane 
    }

    // Allocate project area and assign values
    area[d].allocate(dims(0), dims(1), dims(2));               // the whole plane

    // Following is for a whole part of plane
    Vec3i offset(0); offset(d) = 1;
    Vec3i  start = left_view_slice_start[d]+offset;
    Vec3i    end = left_view_slice_end  [d]+offset;
    int     size = (end-start).product();
    double *data = new double[size];
    grid.area(d).pack(start, end, data);
    for(int i=0; i<size; i++) area[d].data()[i] = data[i];
    delete [] data;

    // Default ranks with periodic setup
    rank_left_vbc(d) = mympi.left(d);
    rank_left_pbc(d) = mympi.left(d);
    rank_left_sbc(d) = mympi.left(d);

    rank_right_vbc(d) = mympi.right(d);
    rank_right_pbc(d) = mympi.right(d);
    rank_right_sbc(d) = mympi.right(d);

  } // loop over 3 directions

  // Allocation for the mean profiles
  if(controller.turbulent_inflow_generation()){

    // Local mean profiles
    umean.allocate(grid.n(2));
    vmean.allocate(grid.n(2));
    wmean.allocate(grid.n(2));

    // Global mean profiles
    gumean.allocate(grid.gn(2));
    gvmean.allocate(grid.gn(2));
    gwmean.allocate(grid.gn(2));

  }

  // Assign values for assistant velocity slices 
  if((controller.get_case() == CASE_URBAN) || (controller.get_case() == CASE_TBL)){
    for(int k=left_view_slice_start[0](2); k<left_view_slice_end[0](2); k++){   // whole domain that involves negative height
      double z=grid.x(2, k)+grid.dx(2, k)/2.-grid.lbox(2);    // get absolute height of stream-wise component
      double U;
      if(z>controller.z0()){
        U = controller.utau()/0.41*log(z/controller.z0());    // log-law profile
      }
      else if (z>0){
        U = z/controller.z0();
      }
      else{
        U = 0;
      }
      for(int j=left_view_slice_start[0](1); j<left_view_slice_end[0](1); j++){
        left_slice[0][0](0,j,k) = U * controller.uref();      // scale by reference velocity 
      }
    }

    // Copy vertical velocity profile for a turbulent inflow generation
    if(controller.turbulent_inflow_generation()){

      // Local mean profile
      {
        int n = 0;
        for(int k=grid.start(2); k<grid.end(2); k++){             // slightly smaller 
          double z=grid.x(2, k)+grid.dx(2, k)/2.-grid.lbox(2);    // get absolute height of stream-wise component
          double U;                                               // THIS IS IDENTICAL TO THE ABOVE ONE
          if(z>controller.z0()){
            U = controller.utau()/0.41*log(z/controller.z0());    // log-law profile
          }
          else if (z>0){
            U = z/controller.z0();
          }
          else{
            U = 0;
          }
          umean(n) = U * controller.uref();
          n++;
        }
      }

      // Global mean profile
      {
        int n = 0;
        for(int k=GC; k<GC+grid.gn(2); k++){         
          double z=(grid.xg(2, k)+grid.xg(2, k+1))/2.-grid.lbox(2);  // get absolute height of stream-wise component
          double U;                                                  // THIS IS IDENTICAL TO THE ABOVE ONE
          if(z>controller.z0()){
            U = controller.utau()/0.41*log(z/controller.z0());       // low-law profile
          }
          else if (z>0){
            U = z/controller.z0();
          }
          else{
            U = 0;
          }
          gumean(n) = U * controller.uref();
          n++;
        }
      }

    } // end if generation

    // Setup specific ranks 
    if(mympi.lboundary(0)){
      rank_left_vbc(0) = -1;
      rank_left_pbc(0) = -1;
      rank_left_sbc(0) = -1;
    }
    if(mympi.lboundary(2)){
      rank_left_vbc(2) = -1;
      rank_left_pbc(2) = -1;
      rank_left_sbc(2) = -1;
    }
    if(mympi.rboundary(0)){
      rank_right_vbc(0) = -1;
      rank_right_pbc(0) = -1;
      rank_right_sbc(0) = -1;
    }
    if(mympi.rboundary(2)){
      rank_right_vbc(2) = -1;
      rank_right_pbc(2) = -1;
      rank_right_sbc(2) = -1;
    }

  } // end if urban case or TBL case

  if(controller.get_case() == CASE_CYLINDER){
    for(int k=left_view_slice_start[0](2); k<left_view_slice_end[0](2); k++){
      for(int j=left_view_slice_start[0](1); j<left_view_slice_end[0](1); j++){
        left_slice[0][0](0,j,k)=controller.uref();
      }
    }

    // Local mean profile
    if(controller.turbulent_inflow_generation())
    {
      int n = 0;
      for(int k=grid.start(2); k<grid.end(2); k++){
        umean(n) = controller.uref();
        n++;
      }
    }

    // Global mean profile
    if(controller.turbulent_inflow_generation())
    {
      int n = 0;
      for(int k=GC; k<GC+grid.gn(2); k++){
        gumean(n) = controller.uref();
        n++;
      }
    }

    // Setup specific ranks 
    if(mympi.lboundary(0)){
      rank_left_vbc(0) = -1;
      rank_left_pbc(0) = -1;
      rank_left_sbc(0) = -1;
    }
    if(mympi.lboundary(2)){
      rank_left_vbc(2) = -1;
      rank_left_pbc(2) = -1;
      rank_left_sbc(2) = -1;
    }
    if(mympi.rboundary(0)){
      rank_right_vbc(0) = -1;
      rank_right_pbc(0) = -1;
      rank_right_sbc(0) = -1;
    }
    if(mympi.rboundary(2)){
      rank_right_vbc(2) = -1;
      rank_right_pbc(2) = -1;
      rank_right_sbc(2) = -1;
    }

  } // end if cylinder case

  if(controller.get_case() == CASE_OPENCHANNEL){
    
    // Setup specific ranks 
    if(mympi.lboundary(2)){
      rank_left_vbc(2) = -1;
      rank_left_pbc(2) = -1;
      rank_left_sbc(2) = -1;
    }
    if(mympi.rboundary(2)){
      rank_right_vbc(2) = -1;
      rank_right_pbc(2) = -1;
      rank_right_sbc(2) = -1;
    }

  } // end if open channel case

  // Initialise turbulent inflow generator
  if(controller.turbulent_inflow_generation()){
    init_turbulent_inflow_generation(grid);
  }
}

// Apply velocity BC based on a specific case
void BC::apply_vbc(VectorFieldd &v,              const MyMPI  &mympi,
                   const Controller &controller, const Grid   &grid ) {

  // Exchange halo cells and apply default periodic bc
  for(int d=0; d<3; d++){
    swap(v[d], mympi, rank_left_vbc, rank_right_vbc);
  }

  // Apply BC assistant variables
  Vec3i get_start, get_end, set_start, set_end, offset;

  // Apply BC with specific case
  int icase = controller.get_case();

  // Apply turbulent inflow generation by adding fluctuations.
  // Here only update left_slice[0][0~3].  
  if(controller.turbulent_inflow_generation()){
    if(mympi.lboundary(0)){
      int kk = 0;          // height index offset
      int n  = 0;          // data index offset
      for(int k=grid.start(2); k<grid.end(2); k++){
        for(int j=grid.start(1); j<grid.end(1); j++){
          left_slice[0][0](0,j,k) = umean(kk) + uf.data()[n];
          left_slice[0][1](0,j,k) = vmean(kk) + vf.data()[n];
          left_slice[0][2](0,j,k) = wmean(kk) + wf.data()[n];
          n++;
        } // end loop j
        kk++;
      } // end loop k
    }
  }

  if((icase == CASE_URBAN) || (icase == CASE_TBL)){

    // Offset in 1st direction for BC in stream-wise direction
    offset = 0; offset(0) = 1;

    // Inlet
    if(mympi.lboundary(0)){
      get_start = left_view_slice_start[0] + offset;
      get_end   = left_view_slice_end  [0] + offset;
      set_start = left_view_slice_start[0];
      set_end   = left_view_slice_end  [0];

      // Treat 1st component and 2,3 components differently
      Dirichlet_bc(v[0], get_start, get_end, left_slice[0][0].data());
      if(controller.turbulent_inflow_generation()) {
        Dirichlet_bc(v[1], set_start, set_end, left_slice[0][1].data());
        Dirichlet_bc(v[2], set_start, set_end, left_slice[0][2].data());
      }
      else {
        Dirichlet_bc(v[1], set_start, set_end, get_start, get_end, left_slice[0][1].data());
        Dirichlet_bc(v[2], set_start, set_end, get_start, get_end, left_slice[0][2].data());
      }

      // u compoment copy -1 slice
      Neumann_bc(v[0], set_start, set_end, get_start, get_end);

      // Copy -2 slices
      for(int d=0; d<3; d++) {
        Neumann_bc(v[d], set_start-offset, set_end-offset, get_start-offset, get_end-offset);
      }
    } // end inlet bc

    // Outlet: convective
    if(mympi.rboundary(0)){
      get_start = right_view_slice_start[0] - offset;
      get_end   = right_view_slice_end  [0] - offset;
      set_start = right_view_slice_start[0];
      set_end   = right_view_slice_end  [0];
      int size = (set_end - set_start).product();
      double *setdata = new double[size];
      double *getdata = new double[size];
      double *Uc      = new double[size];
      double  dt, dx;
      dt = controller.dt();    
      dx = grid.dx(0, grid.last(0));
      v[0].pack(get_start, get_end, Uc);       // get only stream-wise velocity
      for(int i=0; i<size; i++){
        if(Uc[i] < 0) Uc[i] = 0.0;             // filter less than 0 part 
      }
      for(int e=0; e<3; e++){
        v[e].pack(set_start, set_end, setdata);
        v[e].pack(get_start, get_end, getdata);
        for(int i=0; i<size; i++) {
          setdata[i] = (setdata[i]+Uc[i]*dt/dx*getdata[i])/(1+Uc[i]*dt/dx);
        }
        v[e].fill(set_start, set_end, setdata);
      }
      delete [] setdata;
      delete [] getdata;
      delete [] Uc;
    } // end outlet BC

    // Offset in 3rd direction for BC in vertical direction
    offset = 0; offset(2) = 1;

    // Bottom: non-slip wall
    if(mympi.lboundary(2)){
      get_start = left_view_slice_start[2] + offset;
      get_end   = left_view_slice_end  [2] + offset;
      set_start = left_view_slice_start[2];
      set_end   = left_view_slice_end  [2];
      for(int e=0; e<3; e++){
        if(e == 2){
          Dirichlet_bc(v[e], get_start, get_end, 0.0);
        }
        else{
          Dirichlet_bc(v[e], set_start, set_end, get_start, get_end);
        }
      }

      // w component copy -1 slice
      Neumann_bc(v[2], set_start, set_end, get_start, get_end);

      // Copy -2 slices
      for(int d=0; d<3; d++){
        Neumann_bc(v[d], set_start-offset, set_end-offset, get_start-offset, get_end-offset);
      }

    } // end bottom BC

    // Top: free slip
    if(mympi.rboundary(2)){
      get_start = right_view_slice_start[2] - offset;
      get_end   = right_view_slice_end  [2] - offset;
      set_start = right_view_slice_start[2];
      set_end   = right_view_slice_end  [2];
      for(int e=0; e<3; e++){
        if(e == 2){
          Dirichlet_bc(v[e], set_start, set_end, 0.0);
        }
        else{
          Neumann_bc(v[e], set_start, set_end, get_start, get_end);
        }
      }

      // Copy +2 slices
      for(int d=0; d<3; d++){
        Neumann_bc(v[d], set_start+offset, set_end+offset, get_start+offset, get_end+offset);
      }

    } // end top BC

    // Flux correction
    double flux_in  = 0.0;
    double flux_out = 0.0;
    Vec3d  box = grid.rbox() - grid.lbox(); box(0) = 1.; // flatten along the 1st dimension
    double area = box.product();                         // get outlet area
    offset = 0; offset(0) = 1;

    // Compute flux in at the inlet
    if(mympi.lboundary(0)){
      Vec3i start = left_view_slice_start[0] + offset;
      Vec3i end   = left_view_slice_end  [0] + offset;
      int size  = (end - start).product();
      double *u = new double[size];     // inlet velocity
      double *a = new double[size];     // inlet area
      v[0].pack(start, end, u);         // get inlet velocity
      grid.area(0).pack(start, end, a); // get inlet area
      for(int i=0; i<size; i++) flux_in += u[i]*a[i];
      delete [] u;
      delete [] a;
    }
    mympi.sum(flux_in);

    // Compute flux out at the outlet
    if(mympi.rboundary(0)){
      get_start = right_view_slice_start[0] - offset;
      get_end   = right_view_slice_end  [0] - offset;
      set_start = right_view_slice_start[0];
      set_end   = right_view_slice_end  [0];
      int size  = (set_end - set_start).product();
      double *u = new double[size];             // outlet velocity
      double *a = new double[size];             // outlet area
      v[0].pack(set_start, set_end, u);         // get outlet velocity
      grid.area(0).pack(get_start, get_end, a); // get outlet area
      for(int i=0; i<size; i++) flux_out += u[i]*a[i];
      delete [] u;
      delete [] a;
    }
    mympi.sum(flux_out);

    // Correct flux
    double delta = (flux_out - flux_in)/area;
    if(mympi.rboundary(0)){
      Vec3i start = right_view_slice_start[0];
      Vec3i end   = right_view_slice_end  [0];
      int size = (end - start).product();
      double *u = new double[size];
      v[0].pack(start, end, u);                // get outlet slice data
      for(int i=0; i<size; i++) u[i] -= delta; // correction
      v[0].fill(start, end, u);                // inset back
      delete [] u;
    }

  } // CASE_URBAN or TBL

  if(icase == CASE_CYLINDER){

    // Offset in 1st direction for BC in stream-wise direction
    offset = 0; offset(0) = 1;

    // Inlet
    if(mympi.lboundary(0)){
      get_start = left_view_slice_start[0] + offset;
      get_end   = left_view_slice_end  [0] + offset;
      set_start = left_view_slice_start[0];
      set_end   = left_view_slice_end  [0];

      // Treat 1st component and 2,3 components differently
      Dirichlet_bc(v[0], get_start, get_end, left_slice[0][0].data());
      if(controller.turbulent_inflow_generation()) {
        Dirichlet_bc(v[1], set_start, set_end, left_slice[0][1].data());
        Dirichlet_bc(v[2], set_start, set_end, left_slice[0][2].data());
      }
      else {
        Dirichlet_bc(v[1], set_start, set_end, get_start, get_end, left_slice[0][1].data());
        Dirichlet_bc(v[2], set_start, set_end, get_start, get_end, left_slice[0][2].data());
      }

      // u compoment copy -1 slice
      Neumann_bc(v[0], set_start, set_end, get_start, get_end);

      // Copy -2 slices
      for(int d=0; d<3; d++) {
        Neumann_bc(v[d], set_start-offset, set_end-offset, get_start-offset, get_end-offset);
      }

    } // end inlet bc

    // Outlet: convective
    if(mympi.rboundary(0)){
      get_start = right_view_slice_start[0] - offset;
      get_end   = right_view_slice_end  [0] - offset;
      set_start = right_view_slice_start[0];
      set_end   = right_view_slice_end  [0];
      int size = (set_end - set_start).product();
      double *setdata = new double[size];
      double *getdata = new double[size];
      double *Uc      = new double[size];
      double  dt, dx;
      dt = controller.dt();    
      dx = grid.dx(0, grid.last(0));
      v[0].pack(get_start, get_end, Uc);       // get only stream-wise velocity
      for(int i=0; i<size; i++){
        if(Uc[i] < 0) Uc[i] = 0.0;             // filter less than 0 part 
      }
      for(int e=0; e<3; e++){
        v[e].pack(set_start, set_end, setdata);
        v[e].pack(get_start, get_end, getdata);
        for(int i=0; i<size; i++) {
          setdata[i] = (setdata[i]+Uc[i]*dt/dx*getdata[i])/(1+Uc[i]*dt/dx);
        }
        v[e].fill(set_start, set_end, setdata);
      }
      delete [] setdata;
      delete [] getdata;
      delete [] Uc;
    } // end outlet BC

    // Offset in 3rd direction for BC in vertical direction
    offset = 0; offset(2) = 1;

    // Bottom: free slip
    if(mympi.lboundary(2)){
      get_start = left_view_slice_start[2] + offset;
      get_end   = left_view_slice_end  [2] + offset;
      set_start = left_view_slice_start[2];
      set_end   = left_view_slice_end  [2];
      for(int e=0; e<3; e++){
        if(e == 2){
          Dirichlet_bc(v[e], get_start, get_end, 0.0);
        }
        else{
          Neumann_bc(v[e], set_start, set_end, get_start, get_end);
        }
      }

      // w component copy -1 slice
      Neumann_bc(v[2], set_start, set_end, get_start, get_end);

      // Copy -2 slices
      for(int d=0; d<3; d++){
        Neumann_bc(v[d], set_start-offset, set_end-offset, get_start-offset, get_end-offset);
      }

    } // end bottom BC

    // Top: free slip
    if(mympi.rboundary(2)){
      get_start = right_view_slice_start[2] - offset;
      get_end   = right_view_slice_end  [2] - offset;
      set_start = right_view_slice_start[2];
      set_end   = right_view_slice_end  [2];
      for(int e=0; e<3; e++){
        if(e == 2){
          Dirichlet_bc(v[e], set_start, set_end, 0.0);
        }
        else{
          Neumann_bc(v[e], set_start, set_end, get_start, get_end);
        }
      }

      // Copy +2 slices
      for(int d=0; d<3; d++){
        Neumann_bc(v[d], set_start+offset, set_end+offset, get_start+offset, get_end+offset);
      }

    } // end top BC

    // Flux correction
    double flux_in  = 0.0;
    double flux_out = 0.0;
    Vec3d  box = grid.rbox() - grid.lbox(); box(0) = 1.; // flatten along the 1st dimension
    double area = box.product();                         // get outlet area
    offset = 0; offset(0) = 1;

    // Compute flux in at the inlet
    if(mympi.lboundary(0)){
      Vec3i start = left_view_slice_start[0] + offset;
      Vec3i end   = left_view_slice_end  [0] + offset;
      int size  = (end - start).product();
      double *u = new double[size];     // inlet velocity
      double *a = new double[size];     // inlet area
      v[0].pack(start, end, u);         // get inlet velocity
      grid.area(0).pack(start, end, a); // get inlet area
      for(int i=0; i<size; i++) flux_in += u[i]*a[i];
      delete [] u;
      delete [] a;
    }
    mympi.sum(flux_in);

    // Compute flux out at the outlet
    if(mympi.rboundary(0)){
      get_start = right_view_slice_start[0] - offset;
      get_end   = right_view_slice_end  [0] - offset;
      set_start = right_view_slice_start[0];
      set_end   = right_view_slice_end  [0];
      int size  = (set_end - set_start).product();
      double *u = new double[size];             // outlet velocity
      double *a = new double[size];             // outlet area
      v[0].pack(set_start, set_end, u);         // get outlet velocity
      grid.area(0).pack(get_start, get_end, a); // get outlet area
      for(int i=0; i<size; i++) flux_out += u[i]*a[i];
      delete [] u;
      delete [] a;
    }
    mympi.sum(flux_out);

    // Correct flux
    double delta = (flux_out - flux_in)/area;
    if(mympi.rboundary(0)){
      Vec3i start = right_view_slice_start[0];
      Vec3i end   = right_view_slice_end  [0];
      int size = (end - start).product();
      double *u = new double[size];
      v[0].pack(start, end, u);                // get outlet slice data
      for(int i=0; i<size; i++) u[i] -= delta; // correction
      v[0].fill(start, end, u);                // inset back
      delete [] u;
    }

  } // end if cylinder case vbc
  
  if(icase == CASE_OPENCHANNEL){

    // Offset in 3rd direction for BC in vertical direction
    offset = 0; offset(2) = 1;

    // Bottom: non-slip wall
    if(mympi.lboundary(2)){
      get_start = left_view_slice_start[2] + offset;
      get_end   = left_view_slice_end  [2] + offset;
      set_start = left_view_slice_start[2];
      set_end   = left_view_slice_end  [2];
      for(int e=0; e<3; e++){
        if(e == 2){
          Dirichlet_bc(v[e], get_start, get_end, 0.0);
        }
        else{
          Dirichlet_bc(v[e], set_start, set_end, get_start, get_end);
        }
      }

      // w component copy -1 slice
      Neumann_bc(v[2], set_start, set_end, get_start, get_end);

      // Copy -2 slices
      for(int d=0; d<3; d++){
        Neumann_bc(v[d], set_start-offset, set_end-offset, get_start-offset, get_end-offset);
      }

    } // end bottom BC

    // Top: free slip
    if(mympi.rboundary(2)){
      get_start = right_view_slice_start[2] - offset;
      get_end   = right_view_slice_end  [2] - offset;
      set_start = right_view_slice_start[2];
      set_end   = right_view_slice_end  [2];
      for(int e=0; e<3; e++){
        if(e == 2){
          Dirichlet_bc(v[e], set_start, set_end, 0.0);
        }
        else{
          Neumann_bc(v[e], set_start, set_end, get_start, get_end);
        }
      }

      // Copy +2 slices
      for(int d=0; d<3; d++){
        Neumann_bc(v[d], set_start+offset, set_end+offset, get_start+offset, get_end+offset);
      }

    } // end top BC
  
  } // end if open channel case

}

// Set pressure boundary condition
void BC::apply_pbc(Fieldd &f, const MyMPI &mympi, const Controller &controller) const {

  // Exchange halo cells and apply default periodic boundary condition
  swap(f, mympi, rank_left_pbc, rank_right_pbc);

  // Apply BC assistant variables
  Vec3i get_start, get_end, set_start, set_end, offset;

  // Apply BC with a specific case
  for(int d=0; d<3; d++){
    offset = 0; offset(d) = 1;

    // Left boundary
    get_start = left_view_slice_start[d] + offset;
    get_end   = left_view_slice_end  [d] + offset;
    set_start = left_view_slice_start[d];
    set_end   = left_view_slice_end  [d];
    if(controller.left_pbc(d) == NEUMANN  ) { Neumann_bc  (f, set_start, set_end, get_start, get_end); }
    if(controller.left_pbc(d) == DIRICHLET) { Dirichlet_bc(f, set_start, set_end, get_start, get_end); }

    // Right boundary
    get_start = right_view_slice_start[d] - offset;
    get_end   = right_view_slice_end  [d] - offset;
    set_start = right_view_slice_start[d];
    set_end   = right_view_slice_end  [d];
    if(controller.right_pbc(d) == NEUMANN  ) { Neumann_bc  (f, set_start, set_end, get_start, get_end); }
    if(controller.right_pbc(d) == DIRICHLET) { Dirichlet_bc(f, set_start, set_end, get_start, get_end); }
  }
}

void BC::apply_sbc(Fieldd &f, const MyMPI &mympi, const Controller &controller) const {

  // Exchange halo cells and apply default periodic boundary condition
  swap(f, mympi, rank_left_sbc, rank_right_sbc);

  // Apply BC if needed to overwrite periodic BC
  Vec3i get_start, get_end, set_start, set_end, offset;

  for(int d=0; d<3; d++){
    offset = 0; offset(d) = 1;

    // Left boundary   
    get_start = left_view_slice_start[d] + offset;
    get_end   = left_view_slice_end  [d] + offset;
    set_start = left_view_slice_start[d];
    set_end   = left_view_slice_end  [d];
    if(controller.left_sbc(d) == NEUMANN  ) { Neumann_bc  (f, set_start, set_end, get_start, get_end); }
    if(controller.left_sbc(d) == DIRICHLET) { Dirichlet_bc(f, set_start, set_end, get_start, get_end); }

    // Right boundary
    get_start = right_view_slice_start[d] - offset;
    get_end   = right_view_slice_end  [d] - offset;
    set_start = right_view_slice_start[d];
    set_end   = right_view_slice_end  [d];
    if(controller.right_sbc(d) == NEUMANN  ) { Neumann_bc  (f, set_start, set_end, get_start, get_end); }
    if(controller.right_sbc(d) == DIRICHLET) { Dirichlet_bc(f, set_start, set_end, get_start, get_end); }
  }
}


void BC::apply_vic(VectorFieldd &v, const Controller &controller, const Grid &grid) {
  
  // Conventional settings.
  // Apply initial condition by using BC with specific cases
  int icase = controller.get_case();
  if((icase == CASE_URBAN) || (icase == CASE_TBL)){
    for(int k=0; k<v.ke()+GC; k++){
      for(int j=0; j<v.je()+GC; j++){
        for(int i=0; i<v.ie()+GC; i++){
          v[0](i,j,k) = left_slice[0][0](0,j,k);  // Only populate slice values to its downstream
        }
      }
    }
  } // end URBAN and TBL case

  // Velocity initial condition for open channel case with developed laminar channel flow velocity profile
  if(icase == CASE_OPENCHANNEL){
    double nu = controller.nu();           // viscosity
    double G  = controller.bforce(0);      // get body force in stream-wise direction 
    double h  = grid.rbox(2)-grid.lbox(2); // channel height
    double scale = 1./(G*h*h/(2*nu))*controller.uref(); 
    for(int k=0; k<v.ke()+GC; k++){
      double z  = h - (grid.x(2, k) + grid.dx(2, k)/2. - grid.lbox(2));
      double uz = G*h*h/(2*nu)*(1-z*z/(h*h)) * scale;    
      for(int j=0; j<v.je()+GC; j++){
        for(int i=0; i<v.ie()+GC; i++){
          
          // Add initial fluctutations
          v[0](i,j,k) = uz * distribution(gen)*MAGNITUDE + uz; 
          v[1](i,j,k) = uz * distribution(gen)*MAGNITUDE;
          v[2](i,j,k) = uz * distribution(gen)*MAGNITUDE;
        }
      }
    } 
  }

  if(icase == CASE_CYLINDER){
    for(int k=0; k<v.ke()+GC; k++){
      for(int j=0; j<v.je()+GC; j++){
        for(int i=0; i<v.ie()+GC; i++){
          v[0](i,j,k) = controller.uref();         
        }
      }
    }
  } // end CYLINDER case

  // Other IC types

  // Customised settings.
  for(int it=0; it<controller.vic_number(); it++){
    Vec3d start = controller.vic_coords_start(it);
    Vec3d end   = controller.vic_coords_end(it);
    for(int k=v.ks(); k<v.ke(); k++){
      for(int j=v.js(); j<v.je(); j++){
        for(int i=v.is(); i<v.ie(); i++){
          Vec3d coords = grid.index2coords(3, i, j, k);
          bool flag = true;
          for(int ii=0; ii<3; ii++){
            flag = flag && (coords(ii) >= start(ii));
            flag = flag && (coords(ii) <= end(ii));
          }
          if(flag){
            v.set(i, j, k, controller.vic_value(it));   
          }
        }
      }
    }
  } // end customised settings
}

void BC::read_turbulent_inflow_file(Array1d &var, const Grid &grid, const char *filename){

  // Only do if the file exists.
  ifstream file;
  file.open(filename);
  if(file.is_open()){

    // Read the height and value info.
    const int num = 200;       // a big enough number to hold the data along the height
    Array2d profile(num, 2);
    string line;
    int i = 0;
    while(getline(file, line)){
      istringstream str(line);
      double a, b;
      str >> a >> b;
      profile(i, 0) = a;
      profile(i, 1) = b;
      i++;
    }
    while(i<num){
      profile(i, 0) = profile(i-1, 0)+100;  // increment height with a very large value
      profile(i, 1) = profile(i-1, 1);      // with a constant profile
      i++;
    }

    // Assign values.
    int n = 0; // global
    for(int k=GC; k<grid.gn(2)+GC; k++){
      double h = (grid.xg(2, k)+grid.xg(2, k+1))/2. - grid.lbox(2);
      for(int i=1; i<num; i++){
        if( (profile(i-1,0)<=h) && (h<profile(i,0)) ){
          var(n) = profile(i-1,1);
          break;
        }  // end if and break
      }    // end searching for loop
      n++; // increment the variable
    }      // end height for loop 
  }        // end file open
  file.close();
}

void BC::init_turbulent_inflow_generation(const Grid &grid){

  // Local cell size
  int n1 = grid.n(1);
  int n2 = grid.n(2);

  // Local fluctuation slice 
  uf.allocate(n1, n2); 
  vf.allocate(n1, n2);
  wf.allocate(n1, n2);

  // Change to global perspective
  n1 = grid.gn(1);
  n2 = grid.gn(2);

  // Local grid size
  double dy = (grid.rbox(1) - grid.lbox(1))/grid.gn(1);
  double dz = (grid.rbox(2) - grid.lbox(2))/grid.gn(2);

  // Allocation based on grid size and read from file
  uu.allocate(n2);  read_turbulent_inflow_file(uu, grid, "uu.ccc");
  vv.allocate(n2);  read_turbulent_inflow_file(vv, grid, "vv.ccc");
  ww.allocate(n2);  read_turbulent_inflow_file(ww, grid, "ww.ccc");

  uv.allocate(n2);  read_turbulent_inflow_file(uv, grid, "uv.ccc");
  uw.allocate(n2);  read_turbulent_inflow_file(uw, grid, "uw.ccc");
  vw.allocate(n2);  read_turbulent_inflow_file(vw, grid, "vw.ccc");

  Lx.allocate(n2);  read_turbulent_inflow_file(Lx, grid, "Lx.ccc");
  Ly.allocate(n2);  read_turbulent_inflow_file(Ly, grid, "Ly.ccc");
  Lz.allocate(n2);  read_turbulent_inflow_file(Lz, grid, "Lz.ccc");

  Lt.allocate(n2);
  ny.allocate(n2);
  nz.allocate(n2);

  guf.allocate(n1, n2); 
  gvf.allocate(n1, n2);
  gwf.allocate(n1, n2);

  a11.allocate(n2);
  a21.allocate(n2);
  a22.allocate(n2);
  a31.allocate(n2);
  a32.allocate(n2);
  a33.allocate(n2);

  psiu.allocate    (n1, n2);
  psiv.allocate    (n1, n2);
  psiw.allocate    (n1, n2);
  psiu_old.allocate(n1, n2);
  psiv_old.allocate(n1, n2);
  psiw_old.allocate(n1, n2);

  // Calculation along the height
  int max_ny=0, max_nz=0;
  for(int i=0; i<n2; i++) {
    // Lt
    if(gumean(i) > 1.e-6) { Lt(i) = Lx(i)/gumean(i); }  // Calculate time scale by using Lx
    if(Lt(i)     < 1.e-6) { Lt(i) = 1.e-6; }
    Lt(i) /= TURBINLET_INTERVALS;
    // ny and nz
    ny(i) = Ly(i)/dy;
    nz(i) = Lz(i)/dz;
    if(ny(i) < 1) { ny(i) = 1; }
    if(nz(i) < 1) { nz(i) = 1; }
    // Get the maximum one
    if(max_ny < ny(i)) { max_ny=ny(i); }
    if(max_nz < nz(i)) { max_nz=nz(i); }
  }

  int Ny = max_ny*2;
  int Nz = max_nz*2;
  bfilty.allocate(2*Ny+1, n2);
  bfiltz.allocate(2*Nz+1, n2);
  random_u.allocate(2*Ny+n1, 2*Nz+n2);
  random_v.allocate(2*Ny+n1, 2*Nz+n2);
  random_w.allocate(2*Ny+n1, 2*Nz+n2);

  // Calculate b filter
  for(int k=0; k<n2; k++){

    // filter along y axis
    {
      int N = 2*ny(k);
      double by_denom = 0;
      for(int j=0; j<2*N+1; j++){
        double value = exp(-PI*((double)std::abs(j-N))/((double)ny(k)));
        by_denom += value*value;
        bfilty(j, k) = value;
      }
      by_denom = sqrt(by_denom);
      for(int j=0; j<2*N+1; j++){
        bfilty(j, k) /= by_denom;
      }
    }

    // filter along z axis
    {
      int N = 2*nz(k);
      double bz_denom = 0;
      for(int j=0; j<2*N+1; j++){
        double value = exp(-PI*((double)std::abs(j-N))/((double)nz(k)));
        bz_denom += value*value;
        bfiltz(j, k) = value;
      }
      bz_denom = sqrt(bz_denom);
      for(int j=0; j<2*N+1; j++){
        bfiltz(j, k) /= bz_denom;
      }
    }

  } // end loop b filter

  // Calculate Lund matrix
  for(int k=0; k<n2; k++){

    // First column
    if(uu(k)>1.e-6){
      a11(k) = sqrt(uu(k));
      a21(k) = uv(k)/a11(k);
      a31(k) = uw(k)/a11(k);
    }
    else{
      a11(k) = 1.e-8;
      a21(k) = 1.e-8;
      a31(k) = 1.e-8;
    }

    // Second column
    a22(k) = vv(k) - a21(k)*a21(k);
    if(a22(k)>1.e-6){
      a22(k) = sqrt(a22(k));
      a32(k) = (vw(k) - a21(k)*a31(k))/a22(k);
    }
    else{
      a22(k) = 1.e-8;
      a32(k) = 1.e-8;
    }

    // Third column
    a33(k) = ww(k) - a31(k)*a31(k) - a32(k)*a32(k);
    if(a33(k) > 1.e-6)
      a33(k) = sqrt(a33(k));
    else
      a33(k) = 1.e-8;

  } // end loop Lund matrix

  // Initialise Psi slice
  random_gen();
  psiu_old = psiu;
  psiv_old = psiv;
  psiw_old = psiw;
}

void BC::random_gen(){

  // Get random slice 
  int size = random_u.size();
  for(int i=0; i<size; i++) { random_u.data()[i] = distribution(gen); }
  for(int i=0; i<size; i++) { random_v.data()[i] = distribution(gen); }
  for(int i=0; i<size; i++) { random_w.data()[i] = distribution(gen); }

  // Calculate psi plane
  int n1 = guf.Nx();    // horizontal dimension
  int n2 = guf.Ny();    // vertical dimension
  for(int k=0; k<n2; k++){
    int Ny = ny(k)*2;   // filter length in horizontal direction
    int Nz = nz(k)*2;   // filter length in vertical direction
    for(int j=0; j<n1; j++){

      // Convolution
      double pu=0, pv=0, pw=0;
      for(int kk=0; kk<2*Nz+1; kk++){
        for(int jj=0; jj<2*Ny+1; jj++){
          double bij = bfilty(jj, k)*bfiltz(kk, k);
          pu += random_u(jj+j, kk+k)*bij;
          pv += random_v(jj+j, kk+k)*bij;
          pw += random_w(jj+j, kk+k)*bij;
        }
      } // end convolution
      psiu(j, k) = pu;
      psiv(j, k) = pv;
      psiw(j, k) = pw;
    }
  }
}

void BC::fluctuation_gen(const MyMPI &mympi, const Controller &controller, const Grid &grid){

  // Retrieve data. Note global dimension
  double dt = controller.dt();
  int n1    = grid.gn(1);
  int n2    = grid.gn(2);

  // Update Psi slice
  random_gen();

  // Get new Psi slice.
  for(int k=0; k<n2; k++){
    double c1 = exp(-PI*dt/(2.*Lt(k)));
    double c2 = 1-exp(-PI*dt/Lt(k));
    if(c2 > 1.e-6) { c2 = sqrt(c2); }
    else           { c2 = 1.e-8;    }

    for(int j=0; j<n1; j++){
      psiu(j, k) = psiu_old(j, k)*c1 + psiu(j, k)*c2;
      psiv(j, k) = psiv_old(j, k)*c1 + psiv(j, k)*c2;
      psiw(j, k) = psiw_old(j, k)*c2 + psiw(j, k)*c2;
    }
  }

  // Update old Psi slice
  psiu_old = psiu;
  psiv_old = psiv;
  psiw_old = psiw;

  // Get fluctuation slice
  double utotal=0, vtotal=0, wtotal=0;
  int ntotal = n1*n2;
  for(int k=0; k<n2; k++){
    for(int j=0; j<n1; j++){
      double u  = psiu(j, k);
      double v  = psiv(j, k);
      double w  = psiw(j, k);
      guf(j, k) = a11(k)*u;
      gvf(j, k) = a21(k)*u+a22(k)*v;
      gwf(j, k) = a31(k)*u+a32(k)*v+a33(k)*w;
      utotal   += guf(j, k);
      vtotal   += gvf(j, k);
      wtotal   += gwf(j, k);
    }
  }

  // Get delta
  double delu = utotal/ntotal;
  double delv = vtotal/ntotal;
  double delw = wtotal/ntotal;

  // Substract the delta
  for(int k=0; k<n2; k++){
    for(int j=0; j<n1; j++){
      guf(j, k) -= delu;
      gvf(j, k) -= delv;
      gwf(j, k) -= delw;
    }
  }

  // Broadcast the result from process 0
  mympi.bcast(guf.data(), guf.size(), 0);
  mympi.bcast(gvf.data(), gvf.size(), 0);
  mympi.bcast(gwf.data(), gwf.size(), 0);

  // Retrieve local data from global slices
  int o1 = grid.gn(1)/mympi.blocks(1)*mympi.coords(1);  // offset
  int o2 = grid.gn(2)/mympi.blocks(2)*mympi.coords(2);
  for(int k=0; k<grid.n(2); k++){
    for(int j=0; j<grid.n(1); j++){
      int oj = o1 + j;
      int ok = o2 + k;
      uf(j, k) = guf(oj, ok);
      vf(j, k) = gvf(oj, ok);
      wf(j, k) = gwf(oj, ok);
    }
  }
}

// Neumann BC: 0 gradient
void Neumann_bc(Fieldd &f,
                const Vec3i &set_start, const Vec3i &set_end,
                const Vec3i &get_start, const Vec3i &get_end){
  int size = (get_end - get_start).product();
  double *setdata = new double[size];
  double *getdata = new double[size];
  f.pack(get_start, get_end, getdata);
  for(int i=0; i<size; i++) setdata[i] = getdata[i];
  f.fill(set_start, set_end, setdata);
  delete [] setdata;
  delete [] getdata;
}

// Dirichlet BC: for a surface in between cell centres with 0 value 
void Dirichlet_bc(Fieldd &f,
                  const Vec3i &set_start, const Vec3i &set_end,
                  const Vec3i &get_start, const Vec3i &get_end){
  int size = (get_end - get_start).product();
  double *setdata = new double[size];
  double *getdata = new double[size];
  f.pack(get_start, get_end, getdata);
  for(int i=0; i<size; i++) setdata[i] = -getdata[i];
  f.fill(set_start, set_end, setdata);
  delete [] setdata;
  delete [] getdata;
}

// Dirichlet BC for a surface in between cell centres with a specified value
void Dirichlet_bc(Fieldd &f,
                  const Vec3i &set_start, const Vec3i &set_end,
                  const Vec3i &get_start, const Vec3i &get_end,
                  const double *data){
  int size = (get_end - get_start).product();
  double *setdata = new double[size];
  double *getdata = new double[size];
  f.pack(get_start, get_end, getdata);
  for(int i=0; i<size; i++) setdata[i] = 2*data[i]-getdata[i];
  f.fill(set_start, set_end, setdata);
  delete [] setdata;
  delete [] getdata;
}

// Dirichlet BC: specified value on the slice
void Dirichlet_bc(Fieldd &f, const Vec3i &set_start, const Vec3i &set_end,
                  const double *data){
  f.fill(set_start, set_end, data); 
}

// Dirichlet BC: specify a value on the slice
void Dirichlet_bc(Fieldd &f, const Vec3i &set_start, const Vec3i &set_end,
                  const double value){
  int size = (set_end-set_start).product();
  double *data = new double [size];
  for(int i=0; i<size; i++) data[i] = value;
  f.fill(set_start, set_end, data);
  delete [] data;
}
