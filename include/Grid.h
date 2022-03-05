// Grid.h
// Grid class holds global and local coordinate, grid spacing and enquiry.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef GRID_H
#define GRID_H

#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

#include "Vec.h"
#include "Array.h"
#include "MyMPI.h"
#include "Field.h"
#include "parameter.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ios;
using std::string;

class Grid{

  private:

    // Global 
    Array1d      _xg[3];         // global coordinate
    int          _gn[3];         // global dimension for cells
    double     _lbox[3];         // global logical left box (whole domain) vertex coordinate
    double     _rbox[3];         // global logical right box (whole domain) vertex coordinate
    int _piece_start[3];         // piece view start
    int   _piece_end[3];         // piece view end

    // Local
    int           _n[3];         // local dimension for cells
    int       _start[3];         // local start index
    int         _end[3];         // local end index
    int        _size=0 ;         // local number of cells
    Array1d       _x[3];         // local coordinate
    Array1d      _dx[3];         // cell spacing
    Array1d     _dx2[3];         // cross cell spacing
    Array1d     _dxi[3];         // cell spacing inverse
    Array1d    _dxi2[3];         // cross cell spacing inverse

    // Local fields
    Fieldd    _delta   ;         // LES filter length scale
    Fieldd       _dv   ;         // cell volume
    Fieldd     _area[3];         // project areas
    double   _lrange[3];         // partition range on the logical left
    double   _rrange[3];         // partition range on the logical right

  public:

    // Constructor and deconstructor
    Grid(){}
    ~Grid(){}

    // Set up a global uniform grid in direction d [0,1,2] with grid spacing h.  
    void set_uniform_line(const int d, double neg, double pos, double h);

    // Set up a global stretched grid in direction d [0,1,2].
    // The grid start with a uniform region and then the spacing is
    // exponentially increased with a factor.
    void set_stretch_line(const int d, double uneg, double upos, double pos, double h, double factor);

    // Initialise the grid after all 3 global grids have been set.
    void init(const MyMPI &);

    // Get size (number of cells) in the current partition.
    int size() const {return _size;}

    // Get values either with index tuple (i,j,k) or direction d. Return a vector if no argument given.
    double dv          (const int i, const int j, const int k) const {return      _dv(i,j,k);  }
    double delta       (const int i, const int j, const int k) const {return   _delta(i,j,k);  }
    double lrange      (const int d)                           const {return      _lrange[d];  }
    double rrange      (const int d)                           const {return      _rrange[d];  }
    int    start       (const int d)                           const {return       _start[d];  }
    int    end         (const int d)                           const {return         _end[d];  }
    int    piece_start (const int d)                           const {return _piece_start[d];  }
    int    piece_end   (const int d)                           const {return   _piece_end[d];  }

    double xg          (const int d, const int i)              const {return   _xg[d](i); }
    double x           (const int d, const int i)              const {return    _x[d](i); }
    double dx          (const int d, const int i)              const {return   _dx[d](i); }
    double dx2         (const int d, const int i)              const {return  _dx2[d](i); }
    double dxi         (const int d, const int i)              const {return  _dxi[d](i); }
    double dxi2        (const int d, const int i)              const {return _dxi2[d](i); }

    int    gn          (const int d)                           const {return   _gn[d]; }
    int    n           (const int d)                           const {return    _n[d]; }
    double lbox        (const int d)                           const {return _lbox[d]; }
    double rbox        (const int d)                           const {return _rbox[d]; }

    Vec3d  lrange      ()                                      const {return Vec3d(_lrange[0], _lrange[1], _lrange[2]);}
    Vec3d  rrange      ()                                      const {return Vec3d(_rrange[0], _rrange[1], _rrange[2]);}

    Vec3i  start       ()                                      const {return Vec3i(_start[0], _start[1], _start[2]);}
    Vec3i  end         ()                                      const {return Vec3i(  _end[0],   _end[1],   _end[2]);}

    Vec3i  gn          ()                                      const {return Vec3i(_gn[0], _gn[1], _gn[2]);}
    Vec3i  n           ()                                      const {return Vec3i( _n[0],  _n[1],  _n[2]);}

    Vec3d  lbox        ()                                      const {return Vec3d(_lbox[0], _lbox[1], _lbox[2]);}
    Vec3d  rbox        ()                                      const {return Vec3d(_rbox[0], _rbox[1], _rbox[2]);}

    Vec3i  piece_start ()                                      const {return Vec3i(_piece_start[0], _piece_start[1], _piece_start[2]);}
    Vec3i  piece_end   ()                                      const {return Vec3i(  _piece_end[0],   _piece_end[1],   _piece_end[2]);}

    // First and last cell index
    int    first       (const int d)                           const {return GC;       }
    int    last        (const int d)                           const {return GC+_n[d]; }

    // Enquiry the project area with a specified normal direction and index
    double area(const int d, const int i, const int j, const int k) const {return _area[d](i,j,k);  }

    const Fieldd area(const int d) const {return _area[d];}

    // Converter
    // Index to coordinate with a specified direction d.
    // A <0 or >=3 of direction d means cell centre for scalars.
    Vec3d index2coords(const int d, const int i, const int j, const int k) const;

    // Coordinate to the index of cell centre.
    // Need to test if a point is inside the partition first.
    Vec3i coords2index(const Vec3d &coords) const;

    // Enquiry if a point is inside the current partition.
    bool inside(const Vec3d &, const MyMPI &) const;
};

#endif // GRID_H
