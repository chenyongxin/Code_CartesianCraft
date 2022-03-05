// SGS.cc
// Sub-grid scale model for LES.
// Code_CartesianCraft (c) Yongxin Chen

#include <SGS.h>

namespace SGS{

  void compute_SGS(Fieldd &ed, const VectorFieldd &v, const VectorFieldd *S, const Grid &grid, const Controller &controller){
    // TODO: Add options for different SGS models
    // if( controller.SGS() == VREMAN ) Vreman(...);
    Vreman(ed, S, grid);
  }

  void Vreman(Fieldd &ed, const VectorFieldd *S, const Grid &grid){

    double a[3][3], b[3][3];

    for(int k=ed.ks(); k<ed.ke(); k++){
      for(int j=ed.js(); j<ed.je(); j++){
        for(int i=ed.is(); i<ed.ie(); i++){

          // Get transposed velocity gradient tensor.
          for(int ii=0; ii<3; ii++){
            for(int jj=0; jj<3; jj++){
              a[jj][ii] = S[ii][jj](i,j,k);
            }
          }

          // Grid spacing square
          double dx2 = grid.dx(0, i) * grid.dx(0, i);
          double dy2 = grid.dx(1, j) * grid.dx(1, j);
          double dz2 = grid.dx(2, k) * grid.dx(2, k);

          // Compute beta
          for(int ii=0; ii<3; ii++){
            for(int jj=0; jj<3; jj++){
              b[ii][jj] =
                dx2 * a[0][ii] * a[0][jj] + 
                dy2 * a[1][ii] * a[1][jj] +
                dz2 * a[2][ii] * a[2][jj];
            }
          }

          // Compute Bb and a_{ij} a_{ij}
          double bb =
            b[0][0]*b[1][1] - b[0][1]*b[0][1] +
            b[0][0]*b[2][2] - b[0][2]*b[0][2] +
            b[1][1]*b[2][2] - b[1][2]*b[1][2];

          double aa = 0;
          for(int ii=0; ii<3; ii++){
            for(int jj=0; jj<3; jj++){
              aa += a[ii][jj]*a[ii][jj];
            }
          }

          // Compute SGS eddy viscosity
          double bb_aa = bb/aa;
          ed(i,j,k) = 0;
          if(bb/aa > 1.e-8){
            ed(i,j,k) = C_VREMAN * sqrt(bb_aa);
          }

        }
      }
    }
  }

}
