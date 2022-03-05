// VectorField.h
// Container and operator for a 3D vector field.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include <iostream>
#include <fstream>

#include "Vec.h"
#include "Field.h"
#include "Array.h"
#include "MyMPI.h"
#include "parameter.h"

using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;

template<class T>
class VectorField{

  private:

    // Set views for interior domain
    void set_view(){
      view_start =    GC;
      view_end   = _n+GC;
      _size      = 3*_n.product();
    }

    // Data and info
    Field<T>   x,y,z;               // 3 components
    Vec3i         _n;               // shape for interior domain
    int     _size =0;               // size of number of interior cells
    Vec3i view_start;               // view for interior domain
    Vec3i view_end  ;

  public:

    // Constructors and deconstructor
    VectorField(){}

    VectorField(const Vec3i &n):
      x(n), y(n), z(n), _n(n) {
        set_view();
      }

    VectorField(const int nx, const int ny, const int nz):
      x(nx, ny, nz), y(nx, ny, nz), z(nx, ny, nz),
      _n(nx, ny, nz) {
        set_view();
      }

    VectorField(const VectorField<T> &b):
      x(b.x), y(b.y), z(b.z), _n(b._n) {
        set_view();
      }

    ~VectorField(){}

    // Allocation for an empty vector field
    void allocate(const int nx, const int ny, const int nz){
      _n = Vec3i(nx, ny, nz);
      x.allocate(nx, ny, nz);
      y.allocate(nx, ny, nz);
      z.allocate(nx, ny, nz);
      set_view();
    }

    void allocate(const Vec3i &n){
      _n = n;
      x.allocate(n);
      y.allocate(n);
      z.allocate(n);
      set_view();
    }

    // Operator overloading
    // = a vector
    VectorField<T>& operator = (const Vec3<T> b){
      x = b(0);
      y = b(1);
      z = b(2);
      return *this;
    }

    // = another field by reshaping and copying its value
    VectorField<T>& operator = (const VectorField<T> &b){
      _n = b._n ;
      x  = b.x  ;
      y  = b.y  ;
      z  = b.z  ;
      set_view();    // change view
      return *this;
    }

    // [] take ith component 
    Field<T>& operator [] (const int d) {
      if(d == 0){ 
        return x;
      }
      else if(d == 1){
        return y;
      }
      else if(d == 2){
        return z;
      }
      else {
        cerr << "Index error in [] when used VectorField." << endl;
        exit(1);
      }
    }

    const Field<T>& operator [] (const int d) const {
      if(d == 0){
        return x;
      }
      else if(d == 1){
        return y;
      }
      else if(d == 2){
        return z;
      }
      else {
        cerr << "Index error in [] when used VectorField." << endl;
        exit(1);
      }
    }

    // () take (d,i,j,k) data
    T& operator () (const int d, const int i, const int j, const int k){
      return (*this)[d](i,j,k);
    } 

    const T& operator () (const int d, const int i, const int j, const int k) const{
      return (*this)[d](i,j,k);
    }

    // Set values of a vector at (i,j,k)
    void set(const int i, const int j, const int k, const Vec3<T> &v){
      for(int d=0; d<3; d++) (*this)[d](i,j,k) = v(d);
    }

    // Get values at (i,j,k) and return a vector
    Vec3<T> get(const int i, const int j, const int k) const {
      Vec3<T> v;
      for(int d=0; d<3; d++) v(d) = (*this)[d](i,j,k);
      return v;
    }

    // Enquiry
    Vec3i start()      const {return view_start;   }  // get start view
    int   start(int d) const {return view_start(d);}  // get start view in direction d
    Vec3i end  ()      const {return view_end;     }  // get end view
    int   end  (int d) const {return view_end(d);  }  // get end view in direction d

    // Get size and shape
    int size() const {return _size;}
    Vec3i  n() const {return _n;   }

    // Start and end range
    int is() const {return view_start(0);}
    int js() const {return view_start(1);}
    int ks() const {return view_start(2);}

    int ie() const {return view_end(0);}
    int je() const {return view_end(1);}
    int ke() const {return view_end(2);}

    void save(ofstream &file) const {
      x.save(file);
      y.save(file);
      z.save(file);
    }

    void read(ifstream &file){
      x.read(file);
      y.read(file);
      z.read(file);
    }
};

// User defined types
typedef VectorField<int>    VectorFieldi;
typedef VectorField<float>  VectorFieldf;
typedef VectorField<double> VectorFieldd;
typedef VectorField<bool>   VectorFieldb;

#endif // VECTORFIELD_H