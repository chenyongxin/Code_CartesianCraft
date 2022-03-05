// Grid.cc
// Grid class holds global and local coordinate, grid spacing and enquiry.
// Code_CartesianCraft (c) Yongxin Chen

#include "Grid.h"

void Grid::set_uniform_line(const int d, double neg, double pos, double h){

  // Sanity check
  if((d<0) || (d>=3)){
    cerr << "Error: Exceeding dimension of grid!" << endl;
    exit(1);
  }
  if(neg >= pos){
    cerr << "Error: Negative position greater than positive position in line " << d << "!" << endl;
    exit(1);
  }
  if(h <= 0){
    cerr << "Error: Non-positive grid spacing in line " << d << "!" << endl;
    exit(1);
  }

  // Assignment
  double len = pos - neg;   // total length
  int      n = ceil(len/h); // num of cell    
  int ntotal = n+1+2*GC;    // num of cell + 1 for num of vertices + 2 * ghost cells for each side
  pos = neg + n*h;          // overwrite positive location
  _gn[d] = n;
  _xg[d].allocate(ntotal);

  // Get the vertex location
  _xg[d](0) = neg-GC*h;
  for(int i=1; i<ntotal; i++) _xg[d](i) = _xg[d](i-1)+h;

  // Get box vertex coordinate
  _lbox[d] = neg;
  _rbox[d] = pos;
}

void Grid::set_stretch_line(const int d, double uneg, double upos, double pos, double h, double factor){

  // Sanity check
  if((d<0) || (d>=3)){
    cerr << "Error: Exceeding dimension of grid!" << endl;
    exit(1);
  }
  if(uneg >= upos){
    cerr << "Error: Negative uniform position greater than positive uniform position in stretched line " << d << "!" << endl;
    exit(1);
  }
  if(upos >= pos){
    cerr << "Error: Positive stretched end less than positive uniform position in stretched line " << d << "!" << endl;
    exit(1);
  }
  if(h <= 0) {
    cerr << "Error: Non-positive minimum grid spacing in stretched line " << d << "!" << endl;
    exit(1);
  }
  if(factor <= 0) {
    cerr << "Error: Non-positive stretching factor in stretched line " << d << "!" << endl;
    exit(1);
  }

  // Compute and overwrite
  // Uniform part
  double ulen = upos - uneg;    // uniform region length
  int n_uniform = ceil(ulen/h); // number of cells in uniform region
  ulen = n_uniform*h;           // overwrite the length in the uniform region
  upos = uneg + n_uniform*h;    // overwrite the positive uniform position

  // Stretched part
  int n_stretch = 0;            // number of cells in stretched region
  double slen = pos - upos;     // stretched part length
  double s = 0.;                // cumulated stretched length
  double dx;                    // spacing
  while(true){
    n_stretch++;
    dx = h*pow(factor, n_stretch);
    s += dx;
    if(s >= slen) {break;}
  }
  pos = upos + slen;                       // overwrite the positive location
  int ntotal = n_uniform+n_stretch+1+2*GC; // total number of vertices

  // Assignment
  _gn[d] = n_uniform + n_stretch;
  _xg[d].allocate(n_uniform+n_stretch+1+2*GC);

  // Get the vertex location
  // Uniform part
  _xg[d](0) = uneg-GC*h;
  for(int i=1; i<=GC+n_uniform; i++){
    _xg[d](i) = _xg[d](i-1)+h;
  }

  // Stretched part
  int n = 1;
  for(int i=GC+n_uniform+1; i<=GC+n_uniform+n_stretch; i++){
    dx = h*pow(factor, n);
    n++;
    _xg[d](i) = _xg[d](i-1)+dx;
  }

  // Positive ghost cells
  for(int i=GC+n_uniform+n_stretch+1; i<ntotal; i++) {
    _xg[d](i) = _xg[d](i-1)+dx;
  }

  // Get box vertex coordinate
  _lbox[d] = uneg;
  _rbox[d] = pos;
}

void Grid::init(const MyMPI &mympi){

  // View and dimension in each partition
  for(int d=0; d<3; d++){
    _n[d]           = _gn[d]/ mympi.blocks(d);       // get division first
    _piece_start[d] = _n[d] * mympi.coords(d);       // start view, this is fixed
    _piece_end[d]   = _n[d] *(mympi.coords(d)+1);    // end view, get it first
    if(mympi.rboundary(d)){                          // change the last partition.
      _n[d] = _gn[d] - (_n[d]* (mympi.blocks(d)-1)); // correct the last partition's dimension
      _piece_end[d] = _piece_start[d] + _n[d];       // end view for the last partition
    }
  }

  // Calculate size
  _size = _n[0] * _n[1] * _n[2];

  // Grid spacing in each partition
  for(int d=0; d<3; d++){

    // View for loop
    _start[d] = GC;
    _end[d]   = GC+_n[d];

    // Allocation for Array. Need to consider ghost cells
    _x[d].allocate(   _n[d]+2*GC);
    _dx[d].allocate(  _n[d]+2*GC);
    _dx2[d].allocate( _n[d]+2*GC);
    _dxi[d].allocate( _n[d]+2*GC);
    _dxi2[d].allocate(_n[d]+2*GC);
  }

  // Retrieve from global data
  for(int d=0; d<3; d++){
    int n = _gn[d] / mympi.blocks(d);        // first Nd-1 blocks with the same num of cells 
    int offset = n * mympi.coords(d);        // offset in cell

    // Assignment
    for(int i=0; i<_n[d]+2*GC; i++){         // loop over local cells
      int gi = offset+i;                     // convert to global index
      _x[d](i)   =  _xg[d](gi);
      _dx[d](i)  =  _xg[d](gi+1) - _xg[d](gi);
      _dxi[d](i) =  1./_dx[d](i);
    }

    for(int i=1; i<_n[d]+2*GC-1; i++){       // one cell smaller  
      int gi = offset+i;
      _dx2[d](i)  = (_xg[d](gi+1) - _xg[d](gi-1))/2.;
      _dxi2[d](i) =  1./_dx2[d](i);
    }

    // Partition range
    _lrange[d] = _x[d](GC);
    _rrange[d] = _x[d](GC+_n[d]); 
  }

  // Allocate and assign for field, where ghost cells are implicitly implemented
  _delta.allocate(_n[0], _n[1], _n[2]);
  _dv.allocate(_n[0], _n[1], _n[2]);
  for(int k=_dv.ks(); k<_dv.ke(); k++){
    for(int j=_dv.js(); j<_dv.je(); j++){
      for(int i=_dv.is(); i<_dv.ie(); i++){
        _dv(i, j, k) = _dx[0](i) * _dx[1](j) * _dx[2](k);
        _delta(i, j, k) = pow(_dv(i,j,k), 1./3.);
      }
    }
  }

  // Allocate and assign project areas
  for(int d=0; d<3; d++){
    _area[d].allocate(_n[0], _n[1], _n[2]);
  }

  // Get area in each normal direction
  for(int d=0; d<3; d++){
    for(int k=_dv.ks(); k<_dv.ke(); k++){
      for(int j=_dv.js(); j<_dv.je(); j++){
        for(int i=_dv.is(); i<_dv.ie(); i++){
          Vec3i ijk(i,j,k);
          _area[d](i,j,k) = _dv(i,j,k)/_dx[d](ijk(d));
        }
      }
    }
  }

  // Write grid
  if(mympi.master()){
    const char *coords[3] = {"x", "y", "z"};
    for(int d=0; d<3; d++){
      ofstream file;
      file.open( (string(GRID_DIR)+string("/")+string(coords[d])).c_str() );
      for(int i=0; i<_gn[d]+1; i++){
        file << _xg[d](GC+i) << endl; 
      }
      file.close();
    }
  }

}

bool Grid::inside(const Vec3d &b, const MyMPI &mympi) const {
  bool a = true;
  for(int d=0; d<3; d++) {
    a = a && (_lrange[d] <= b(d));
    if(mympi.rboundary(d)) {
      a = a && (b(d) <= _rrange[d]);
    }
    else{
      a = a && (b(d) < _rrange[d]);
    }
  }
  return a;
}

Vec3d Grid::index2coords(const int d, const int i, const int j, const int k) const {
  Vec3d coords;
  int ijk[3] = {i, j, k};
  for(int d=0; d<3; d++){
    coords(d) = 0.5*( _x[d](ijk[d]) + _x[d](ijk[d]+1) ); // cell centre coords
  }
  if((0<=d) && (d<3)){
    coords(d) = _x[d](ijk[d]);  // use surface coordinate to overwrite it
  }
  return coords;
}

Vec3i Grid::coords2index(const Vec3d &coords) const {
  Vec3i index;
  for(int d=0; d<3; d++){
    index(d) = _start[d];     // init it with starting index
  }
  for(int d=0; d<3; d++){
    int i = index(d);    // get it first
    while( 0.5*( _x[d](i) + _x[d](i+1) ) < coords(d) ){
      i++;
      if(i >= _end[d]) break;
    }
    index(d) = i-1;      // give it back then
  }
  return index;
}
