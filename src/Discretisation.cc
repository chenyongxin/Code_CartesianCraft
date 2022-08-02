// Discretisation.cc
// Namespace to discretise the governing equations.
// Code_CartesianCraft (c) Yongxin Chen

#include "Discretisation.h"
namespace Discretisation{

  double GammaScheme(double uj, double vm1, double v, double vp1, double vp2);

  void project(VectorFieldd &v, const Fieldd &p, const Grid &grid, const Controller &controller){

    // Get time step.
    double dt = controller.dt();

    // Loop over components and make piecewise update.
    for(int k=p.ks(); k<p.ke(); k++){
      for(int j=p.js(); j<p.je(); j++){
        for(int i=p.is(); i<p.ie(); i++){
          v[0](i,j,k) -= grid.dxi2(0, i)*(p(i, j, k) - p(i-1, j, k))*dt;
          v[1](i,j,k) -= grid.dxi2(1, j)*(p(i, j, k) - p(i, j-1, k))*dt;
          v[2](i,j,k) -= grid.dxi2(2, k)*(p(i, j, k) - p(i, j, k-1))*dt;
        }
      }
    }
  }

  void divergence(Fieldd &f, const VectorFieldd &v, const Grid &grid){
    for(int k=v.ks(); k<v.ke(); k++){
      for(int j=v.js(); j<v.je(); j++){
        for(int i=v.is(); i<v.ie(); i++){
          f(i,j,k) = grid.dv(i,j,k) * ((v[0](i+1, j, k) - v[0](i, j, k)) * grid.dxi(0, i) +
              (v[1](i, j+1, k) - v[1](i, j, k)) * grid.dxi(1, j) +
              (v[2](i, j, k+1) - v[2](i, j, k)) * grid.dxi(2, k) );
        }
      }
    }
  }

  void gradtensor(VectorFieldd *S_stag, VectorFieldd *S_cen, const VectorFieldd &v, const Grid &grid){
    Vec3i o, si, sj; // offset, offset in i and j
    int is = v.is(), js = v.js(), ks = v.ks();
    int ie = v.ie(), je = v.je(), ke = v.ke();

    for(int d1=0; d1<3; d1++){
      for(int d2=0; d2<3; d2++){

        // Shift
        o =0;  o(d1)=1; o(d2)=1;
        si=0; si(d1)=1;
        sj=0; sj(d2)=1;

        if(d1 == d2){    // diagonal component
          for(int k=ks-o(2); k<ke; k++){
            for(int j=js-o(1); j<je; j++){
              for(int i=is-o(0); i<ie; i++){
                Vec3i ijk(i,j,k);
                S_stag[d1][d2](i,j,k) = grid.dxi(d2, ijk(d2))*(v[d1](i+sj(0),j+sj(1),k+sj(2))-
                                                               v[d1](i,      j,      k      ));
              }
            }
          }
          for(int k=ks; k<ke; k++){
            for(int j=js; j<je; j++){
              for(int i=is; i<ie; i++){
                S_cen[d1][d2](i,j,k)  = S_stag[d1][d2](i,j,k);
              }
            }
          }
        }
        else{           // off-diagonal component
          for(int k=ks; k<ke+o(2); k++){
            for(int j=js; j<je+o(1); j++){
              for(int i=is; i<ie+o(0); i++){
                Vec3i ijk(i,j,k);
                S_stag[d1][d2](i,j,k) = grid.dxi2(d2, ijk(d2))*(v[d1](i,      j,      k      )-
                                                                v[d1](i-sj(0),j-sj(1),k-sj(2)));
              }
            }
          }
          for(int k=ks; k<ke; k++){
            for(int j=js; j<je; j++){
              for(int i=is; i<ie; i++){
                S_cen[d1][d2](i,j,k) = 0.25*(S_stag[d1][d2](i,      j,      k     ) +
                                             S_stag[d1][d2](i+si(0),j+si(1),k+si(2))+
                                             S_stag[d1][d2](i+sj(0),j+sj(1),k+sj(2))+
                                             S_stag[d1][d2](i+ o(0),j+ o(1),k+ o(2)));
              }
            }
          }
        } // end if
      } // end for d2
    } // end for d1
  }

  void conv_diff(VectorFieldd &a, const VectorFieldd &u, const Grid &grid, const double nu, const Fieldd &ed, const VectorFieldd *uij){
    a = Vec3d(0.0);          // reinit the result vector field
    Vec3i o, si, sj;         // offset, offset in i and j
    Fieldd effvs(grid.n());  // temporary effective viscosity
    Fieldd flux(grid.n());   // temporary flux
    int is = u.is(), js = u.js(), ks = u.ks();
    int ie = u.ie(), je = u.je(), ke = u.ke();

    // Convective term
    for(int d1=0; d1<3; d1++){
      for(int d2=0; d2<3; d2++){

        // Shift
        o =0;  o(d1)=1; o(d2)=1;
        si=0; si(d1)=1;
        sj=0; sj(d2)=1;

        if(d1 == d2){ // diagonal

          // Compute effective viscosity
          for(int k=ks-o(2); k<ke; k++){
            for(int j=js-o(1); j<je; j++){
              for(int i=is-o(0); i<ie; i++){
                effvs(i,j,k) = ed(i,j,k)+nu;
              }
            }
          }

          // Compute diffusion term
          for(int k=ks; k<ke; k++){
            for(int j=js; j<je; j++){
              for(int i=is; i<ie; i++){
                Vec3i ijk(i,j,k);
                double p0  = uij[d1][d2](i,      j,      k      )+uij[d2][d1](i,      j,      k      );
                double pm1 = uij[d1][d2](i-sj(0),j-sj(1),k-sj(2))+uij[d2][d1](i-sj(0),j-sj(1),k-sj(2));
                a[d1](i,j,k) += grid.dxi2(d2, ijk(d2))*(effvs(i,j,k)*p0-effvs(i-sj(0),j-sj(1),k-sj(2))*pm1);
              }
            }
          }
        }
        else{        // off-diagonal

          // Compute effective viscosity by linear interpolation at edge centre
          for(int k=ks; k<ke+o(2); k++){
            for(int j=js; j<je+o(1); j++){
              for(int i=is; i<ie+o(0); i++){
                effvs(i,j,k) = 0.25*(ed(i,      j,      k      ) +
                                     ed(i-si(0),j-si(1),k-si(2)) +
                                     ed(i-sj(0),j-sj(1),k-sj(2)) +
                                     ed(i- o(0),j- o(1),k- o(2)));
              }
            }
          }

          // Compute convective term
          for(int k=ks; k<ke; k++){
            for(int j=js; j<je; j++){
              for(int i=is; i<ie; i++){
                Vec3i ijk(i,j,k);
                double p0  = uij[d1][d2](i,      j,      k      )+uij[d2][d1](i,      j,      k      );
                double pp1 = uij[d1][d2](i+sj(0),j+sj(1),k+sj(2))+uij[d2][d1](i+sj(0),j+sj(1),k+sj(2));
                a[d1](i,j,k) += grid.dxi(d2, ijk(d2))*(effvs(i+sj(0),j+sj(1),k+sj(2))*pp1 -
                                                       effvs(i,      j,      k      )*p0) ;
              }
            }
          }
        } // end diagonal or off-diagonal
      } // end d2 loop
    } // end d1 loop

    // Compute convective term
    for(int d1=0; d1<3; d1++){
      for(int d2=0; d2<3; d2++){

        // Shift
        o =0;  o(d1)=1; o(d2)=1;
        si=0; si(d1)=1;
        sj=0; sj(d2)=1;

        for(int k=ks; k<ke+o(2); k++){
          for(int j=js; j<je+o(1); j++){
            for(int i=is; i<ie+o(0); i++){
              double uj = 0.5*(u[d2](i,j,k)+u[d2](i-si(0),j-si(1),k-si(2)));
              //double ui = 0.5*(u[d1](i,j,k)+u[d1](i-sj(0),j-sj(1),k-sj(2)));  // Central difference scheme
              
              // Bounded Gamma scheme from Jasak 1996 PhD thesis.   
              double ui = GammaScheme(uj, u[d1](i-2*sj(0), j-2*sj(1), k-2*sj(2)),
                                          u[d1](i-  sj(0), j-  sj(1), k-  sj(2)),
                                          u[d1](i,         j,         k        ),
                                          u[d1](i+  sj(0), j+  sj(1), k+  sj(2)));
              
              flux(i,j,k) = ui*uj;
            }
          }
        }
        if(d1 == d2){
          for(int k=ks; k<ke; k++){
            for(int j=js; j<je; j++){
              for(int i=is; i<ie; i++){
                Vec3i ijk(i,j,k);
                a[d1](i,j,k) -= grid.dxi2(d2,ijk(d2))*(flux(i+sj(0),j+sj(1),k+sj(2))-
                                                       flux(i,      j,      k      ));
              }
            }
          }
        }
        else{
          for(int k=ks; k<ke; k++){
            for(int j=js; j<je; j++){
              for(int i=is; i<ie; i++){
                Vec3i ijk(i,j,k);
                a[d1](i,j,k) -= grid.dxi(d2,ijk(d2))*(flux(i+sj(0),j+sj(1),k+sj(2))-
                                                      flux(i,      j,      k      ));
              }
            }
          }
        } 
      } // end d2
    } // end d1
  } // end conv_diff

  void AB(VectorFieldd &v, const VectorFieldd &r, const VectorFieldd &h, const Vec3d &bforce, const double dt){
    for(int d=0; d<3; d++){
      for(int k=v.ks(); k<v.ke(); k++){
        for(int j=v.js(); j<v.je(); j++){
          for(int i=v.is(); i<v.ie(); i++){
            v[d](i,j,k) += (1.5*r[d](i,j,k) - 0.5*h[d](i,j,k) + bforce(d))*dt;
          }
        }
      }
    }
  }

  // Compute drag force D_i = -Cd u_i V, where V=sqrt(u_i u_i)
  void get_dragforce(const VectorFieldd &u, VectorFieldd &D, const Grid &grid,  const double Cd, const double h){

    int is = u.is(), js = u.js(), ks = u.ks();
    int ie = u.ie(), je = u.je(), ke = u.ke();
    
    for(int d=0; d<3; d++){
      for(int k=ks; k<ke; k++){
        
        // Get current height
        double height;
        if(d == 2){
          height = grid.x(2, k+1);
        }
        else{
          height = grid.x(2, k) + 0.5*grid.dx(2, k);
        }

        // Below drag canopy height
        if(height <= h){
          for(int j=js; j<je; j++){
            for(int i=is; i<ie; i++){
              double ui = 0.5*(u[0](i, j, k)+u[0](i+1, j, k));
              double vi = 0.5*(u[1](i, j, k)+u[1](i, j+1, k));
              double wi = 0.5*(u[2](i, j, k)+u[2](i, j, k+1));
              D[d](i,j,k) = -Cd*u[d](i,j,k)*sqrt(ui*ui+vi*vi+wi*wi);
            } // i loop
          } // j loop
        } // height <= h
      } // k loop 
    } // direction loop

  }

  // u = u + other
  void add_vectorfields(VectorFieldd &u, const VectorFieldd &other){
    int is = u.is(), js = u.js(), ks = u.ks();
    int ie = u.ie(), je = u.je(), ke = u.ke();
    for(int d=0; d<3; d++){
      for(int k=ks; k<ke; k++){
        for(int j=js; j<je; j++){
          for(int i=is; i<ie; i++){
            u[d](i,j,k) += other[d](i,j,k);  
          }
        }
      }  
    }
  }
  


  double GammaScheme(double uj, double vm1, double v, double vp1, double vp2){
    double u, c, d, phic, r;   // upper-stream, current cell, down-stream
    double result;
    if(uj >= 0){
      c = v;
      u = vm1;                 // v - 1
      d = vp1;                 // v + 1
    }    
    else{
      c = vp1;
      u = vp2;
      d = v;
    }
    if(d == u){
      phic = (c-u)/(1.e-8);
    }
    else{
      phic = (c-u)/(d-u);
    }
    if( (phic <= 0) || (phic >= 1) ){
      result = c;         // upwind scheme
    } 
    if( (phic >= BETA) && (phic < 1) ){
      result = 0.5*(c+d); // central difference scheme
    }
    if( (phic > 0) && (phic < BETA) ){
      double r = phic/BETA;
      result = (1-0.5*r)*c+0.5*r*d;  // blended
    }
    return result;
  }

} // end namespace Discretisation

