// Field.h
// Container and operator for a 3D field.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef FIELD_H
#define FIELD_H

#include <fstream>
#include "Vec.h"
#include "Array.h"
#include "MyMPI.h"
#include "parameter.h"

using std::ifstream;
using std::ofstream;

// Class template for 3D field
template<class T>
class Field {

  private:

    // Set views for interior domain. Called followed by shape initalised
    void set_view(){
      view_start =    GC;
      view_end   = _n+GC;
      _size      = _n.product();
    }

    // Data and info
    Array::Array3<T> p;                 // data
    Vec3i           _n;                 // shape for interior domain
    int        _size=0;                 // size of number of interior cells

    // View: Indices [start, end) for the part of a 3D field.
    Vec3i view_start;                   // view for interior domain
    Vec3i view_end;

  public:

    // Constructors and deconstructor
    Field(){}

    Field(const Vec3i &n):
      p(n(0)+2*GC, n(1)+2*GC, n(2)+2*GC),
      _n(n) {
        set_view();
      }

    Field(const int nx, const int ny, const int nz):
      p(nx+2*GC, ny+2*GC, nz+2*GC),
      _n(Vec3i(nx, ny, nz)) {
        set_view();
      }

    Field(const Field<T> &b):
      p(b.p),
      _n(b._n) {
        set_view();
      }

    ~Field(){}


    // Allocation for an empty field
    void allocate(const int nx, const int ny, const int nz){
      _n = Vec3i(nx, ny, nz);
      p.allocate(nx+2*GC, ny+2*GC, nz+2*GC);
      set_view();
    }

    void allocate(const Vec3i &n){
      _n = n;
      p.allocate(n(0)+2*GC, n(1)+2*GC, n(2)+2*GC);
      set_view();
    }

    // Operator overloading
    // = a value
    Field<T>& operator = (const T b){
      p = b;         // assign a value to an allocated array and field
      return *this;
    }

    // = another field by reshaping and copying its value
    Field<T>& operator = (const Field<T> &b){
      _n = b._n;
      p  = b.p;
      set_view();    // change view
      return *this;
    }

    // (i,j,k) 
    T& operator () (const int i, const int j, const int k) {
      return p(i,j,k);
    }

    const T& operator () (const int i, const int j, const int k) const {
      return p(i,j,k);
    }

    // Fill data to field with specified view
    void fill(const Vec3i &start, const Vec3i &end, const T *data){
      int n = 0;
      for(int k=start(2); k<end(2); k++){
        for(int j=start(1); j<end(1); j++){
          for(int i=start(0); i<end(0); i++){
            p(i,j,k) = data[n];
            n++;
          }
        }
      }
    }

    // Pack (copy) data from field with specified view
    void pack(const Vec3i &start, const Vec3i &end, T *data) const {
      int n = 0;
      for(int k=start(2); k<end(2); k++){
        for(int j=start(1); j<end(1); j++){
          for(int i=start(0); i<end(0); i++){
            data[n] = p(i,j,k);
            n++;
          }
        }
      }
    }

    Vec3i start()      const {return view_start;   }  // get start view
    int   start(int d) const {return view_start(d);}  // get start view in direction d
    Vec3i end  ()      const {return view_end;     }  // get end view
    int   end  (int d) const {return view_end(d);  }  // get end view in direction d

    // Get max
    T max(const MyMPI &mympi) const {
      T a = -1e10;
      for(int k=start(2); k<end(2); k++){
        for(int j=start(1); j<end(1); j++){
          for(int i=start(0); i<end(0); i++){
            if(a < p(i,j,k)) a = p(i,j,k);
          }
        }
      }
      mympi.max<T>(a);
      return a;
    }

    // Get min
    T min(const MyMPI &mympi) const {
      T a = 1e10;
      for(int k=start(2); k<end(2); k++){
        for(int j=start(1); j<end(1); j++){
          for(int i=start(0); i<end(0); i++){
            if(p(i,j,k) < a) { a = p(i,j,k); }
          }
        }
      }
      mympi.min<T>(a);
      return a;
    }

    // Scale a field
    void scale(double factor){
      for(int k=start(2); k<end(2); k++){
        for(int j=start(1); j<end(1); j++){
          for(int i=start(0); i<end(0); i++){
            p(i,j,k) *= factor;
          }
        }
      }
    }

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

    // Take data
    Array::Array3<T>& data() {return p;}
    const Array::Array3<T>& data() const {return p;}
    
    // Save binary data with binary file object
    void save(ofstream &file) const {
      for(int k=0; k<_n(2)+2*GC; k++){
        for(int j=0; j<_n(1)+2*GC; j++){
          for(int i=0; i<_n(0)+2*GC; i++){
            double x = (double)p(i, j, k);
            file.write((char*)&x, sizeof(double));
          }
        }
      }
    }

    // Read binary data with binary file object
    void read(ifstream &file){
      double x;
      for(int k=0; k<_n(2)+2*GC; k++){
        for(int j=0; j<_n(1)+2*GC; j++){
          for(int i=0; i<_n(0)+2*GC; i++){
            file.read((char*)&x, sizeof(double));
            p(i, j, k) = x;
          }
        }
      }
    }
};

// User defined types
typedef Field<int>    Fieldi;
typedef Field<float>  Fieldf;
typedef Field<double> Fieldd;
typedef Field<bool>   Fieldb;

#endif // FIELD_H