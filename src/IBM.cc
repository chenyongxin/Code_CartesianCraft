// IBM.cc
// Immersed boundary method in 3D.
// Code_CartesianCraft (c) Yongxin Chen

#include "IBM.h"

void IBM::init(const Grid &grid, const Controller &controller, const MyMPI &mympi){

  // Allocate fields if IBM is needed
  if(controller.get_ibmtype() != 0){
    for(int d=0; d<4; d++){
      _distance[d].allocate(grid.n());
      _normal[d].allocate(grid.n());
      _mask[d].allocate(grid.n());
    }
  }

  // Read geometry input file
  if(controller.get_ibmtype() == CUBOIDS){
    ifstream file("cuboids.ccc");
    if(!file.is_open()){
      if(mympi.rank()==0){ cerr << "Input file ``cuboids.ccc'' does not exist." << endl;}
      exit(1);
    }

    // Get number of cuboids and allocate 2D arrays
    string line;
    getline(file, line);
    istringstream(line) >> ncuboids;
    xc.allocate(ncuboids, 3);
    dh.allocate(ncuboids, 3);

    // Read geometry info and assign values
    double value;
    string values;
    for(int i=0; i<ncuboids; i++){
      getline(file, line);
      istringstream values(line);
      for(int j=0; j<3; j++){
        values >> value;
        xc(i,j) = value;
      }
      for(int j=0; j<3; j++){
        values >> value;
        dh(i,j) = value;
      }
    }

    // Compute wall normal and distance
    distance_normal_aligned_cuboids(ncuboids, grid, xc, dh, _distance, _normal);
  }
  else if(controller.get_ibmtype() == URBAN){
    if(mympi.rank()==0) cerr << "Urban roughness module is under development." << endl;
    exit(1);
  }

  // Compute mask arrays to find the first grid point
  if(controller.get_ibmtype() != 0){
    for(int d=0; d<4; d++){
      for(int k=grid.start(2); k<grid.end(2); k++){
        for(int j=grid.start(1); j<grid.end(1); j++){
          for(int i=grid.start(0); i<grid.end(0); i++){
            if(((_distance[d](i-1,j,k)<0) || (_distance[d](i+1,j,k)<0)  ||
                (_distance[d](i,j-1,k)<0) || (_distance[d](i,j+1,k)<0)  ||
                (_distance[d](i,j,k-1)<0) || (_distance[d](i,j,k+1)<0)) &&
                (_distance[d](i,j,k)>=0)) _mask[d](i,j,k) = true;
          }
        }
      }
    }
  }
}


void IBM::ibc (VectorFieldd &v) const {

  // Only correct 3 velocity components
  for(int d=0; d<3; d++){
    for(int k=v.ks(); k<v.ke(); k++){
      for(int j=v.js(); j<v.je(); j++){
        for(int i=v.is(); i<v.ie(); i++){
          if(distance(d,i,j,k)<=0){   // inside the body
            v[d](i,j,k) = 0.;
            continue;
          }
          if(mask(d,i,j,k)){
            Vec3i o;                  // offset
            for(int ii=0; ii<3; ii++) o(ii) = round(normal(d,i,j,k)(ii));
            double d1 = distance(d, i,      j,      k     );
            double d2 = distance(d, i+o(0), j+o(1), k+o(2));
            v[d](i,j,k) = v[d](i+o(0), j+o(1), k+o(2)) - v[d](i+o(0), j+o(1), k+o(2))/d2*(d2-d1);
          }
        }  
      }
    }
  }
}


// Get distance from a point `pt` to a plane with centre coords `xc` and normal vector `normal`
double distance2plane(const Vec3d &pt, const Vec3d &xc, const Vec3d &normal){
  Vec3d v = pt - xc;
  return v.dot(normal);
}

void distance_normal_aligned_cuboid(const Vec3d &pt, const Vec3d &xc, const Vec3d &dh,
    double &distance, Vec3d &normal){
  double dist;             // temporary distance
  distance = -1.e5;        // init distance
  Vec3d o, cen, n;         // offset, centre coordinate and normal vector
  for(int d=0; d<3; d++){
    o = 0.; o(d) = 1.; 

    // left plane: dh is side length, should divide 2 for half length.
    cen = xc - o*dh/2.;
    n = 0.; n(d) = -1.;
    dist = distance2plane(pt, cen, n);
    if(dist > distance){
      distance = dist;
      normal = n;
    }

    // right plane
    cen = xc + o*dh/2.;
    n = 0.; n(d) = 1.;
    dist = distance2plane(pt, cen, n);
    if(dist > distance){
      distance = dist;
      normal = n;
    }
  }
}

void distance_normal_aligned_cuboids(int n, const Grid &grid, const Array2d &xc,
    const Array2d &dh, Fieldd *distance, VectorFieldd *normal){
  Vec3d pt, norm, newnorm, side, centre;

  // 4 components
  for(int d=0; d<4; d++){ 

    // Loop each point. Get 1 more points in each side of the domain.
    for(int k=grid.start(2)-1; k<grid.end(2)+1; k++){
      for(int j=grid.start(1)-1; j<grid.end(1)+1; j++){
        for(int i=grid.start(0)-1; i<grid.end(0)+1; i++){

          // Get current point.
          pt = grid.index2coords(d, i, j, k);
          double dist=1.e5, newdist;    // reset in every point

          // Loop over all cuboids. 
          for(int ii=0; ii<n; ii++){    

            // Get data: centre and side length vector of each cuboid.
            for(int jj=0; jj<3; jj++){  
              centre(jj) = xc(ii, jj);
              side(jj) = dh(ii, jj);
            }

            // Compute distance and normal direction, and assign data.
            // Find the closest distance and normal direction.
            distance_normal_aligned_cuboid(pt, centre, side, newdist, newnorm);
            
            if(newdist < dist){
              dist = newdist;
              norm = newnorm;
            }
          } // end loop all cuboids

          // Transfer data back
          distance[d](i,j,k) = dist;
          normal[d].set(i, j, k, norm);
        }
      }
    } // end each point 
  } // end 4 components
}
