// Array.h
// Class template for array up to three-dimension.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>

namespace Array{

  const int unallocated = 0;    // status

  inline void ArrayExit(const char* x){
    std::cerr << std::endl << "Error: " << x << "." << std::endl; 
    exit(1);
  }

  // 1D array
  template<class T>
    class Array1{

      private:
        T  *v;                    // data
        int _size;                // total size
        int status;               // allocation status

      public:

        // Constructors and deconstructor
        Array1(): v(NULL), _size(0), status(unallocated) {}

        Array1(const int n): v(NULL), _size(0), status(unallocated){
          allocate(n);
        }

        Array1(const Array1<T> &b): v(NULL), _size(0), status(unallocated){
          allocate(b._size);
          for(int i=0; i<b._size; i++){ 
            v[i] = b.v[i]; 
          }
        }

        ~Array1(){
          if(allocated()){ 
            deallocate(); 
          }
        }

        // Allocation related
        void allocate(const int n){
          if(allocated()){ 
            ArrayExit("Try to allocate an allocated Array1"); 
          }
          if(n<=0){ 
            ArrayExit("Try to allocate non-positive size to Array1"); 
          }
          _size = n;
          v = new T[_size];
          for(int i=0; i<n; i++){ 
            v[i] = 0; 
          }
          status = 1;        	
        }

        void deallocate(){
          if(!allocated()){
            ArrayExit("Try to deallocate a unallocated Array 1");
          }
          delete [] v;
          _size = 0;
          status = unallocated;
        }

        bool allocated() {return status && (_size != 0);}  

        // Enquiry 
        int Nx(){return _size;}

        // Operator overwriting
        Array1<T>& operator = (const Array1<T> &b){
          if(allocated()){
            if(_size != b._size){
              deallocate();
              allocate(b._size);
            }
          }
          else{
            allocate(b._size);
          }
          for(int i=0; i<b._size; i++) {
            v[i] = b.v[i]; 
          }
          return *this;
        }

        Array1<T>& operator = (const T b){
          if(!allocated()){
            ArrayExit("Try to assign a value to a unallocated Array1");
          }
          for(int i=0; i<_size; i++) { 
            v[i] = b; 
          }
          return *this;
        }

        T& operator () (const int i) {
          if((i<0) || (i>=_size)) {
            ArrayExit("Array1 index out of bound");
          }
          return v[i];
        }

        const T& operator () (const int i) const {
          if((i<0) || (i>=_size)) {
            ArrayExit("Array1 index out of bound");
          }
          return v[i];
        }

        // Take data
        T* data() {return v;}
        const T* data() const {return v;}

        // Get size
        int size() const {return _size;}
    };

  // 2D array
  template<class T> 
    class Array2{

      private:
        T  *v;                // data
        int _size;            // total size
        int status;           // allocation status
        int nx = 0;           // 1st dimension
        int ny = 0;           // 2nd dimension

      public:

        // Constructor and deconstructor
        Array2(): v(NULL), _size(0), status(unallocated) {}

        Array2(const int nx0, const int ny0): v(NULL), _size(0), status(unallocated) {
          allocate(nx0, ny0);
        }

        Array2(const Array2<T> &b): v(NULL), _size(0), status(unallocated) {
          allocate(b.nx, b.ny);
          for(int i=0; i<b._size; i++) {
            v[i] = b.v[i];
          }
        }

        ~Array2(){
          if(allocated()){
            deallocate();
          }
        }

        // Allocation and deallocation
        void allocate(int nx0, int ny0){
          if(allocated()){
            ArrayExit("Try to allocate an allocated Array2");
          }
          if((nx0 <= 0) || (ny0 <= 0)){
            ArrayExit("Try to allocate non-positive size in Array2");
          }
          nx = nx0;
          ny = ny0;
          _size = nx * ny;
          v = new T[_size];
          for(int i=0; i <_size; i++){
            v[i] = 0;
          }
          status = 1;
        }

        void deallocate(){
          if(!allocated()){
            ArrayExit("Try to deallocate a unallocated Array 2");
          }
          delete [] v;
          _size = 0;
          nx = 0;
          ny = 0;
          status = unallocated;
        }

        bool allocated(){return status && (_size != 0);}  

        // Enquiry
        int Nx() const {return nx;}
        int Ny() const {return ny;}

        // Operator overwriting
        Array2<T>& operator = (const Array2<T> &b){
          if(allocated()){
            if((nx != b.nx) || (ny != b.ny)){
              deallocate();
              allocate(b.nx, b.ny);
            }
          }
          else{
            allocate(b.nx, b.ny);
          }
          for(int i=0; i<b._size; i++){
            v[i] = b.v[i];
          }
          return *this;
        }

        Array2<T>& operator = (const T b){
          if(!allocated()){
            ArrayExit("Try to assign a value to a unallocated Array2");
          }
          for(int i=0; i<_size; i++){
            v[i] = b;
          }
          return *this;
        }

        T& operator () (const int i, const int j) {
          if((i<0) || (i>=nx)){
            ArrayExit("Array2 1st index out of bound");
          } 
          if((j<0) || (j>=ny)){
            ArrayExit("Array2 2nd index out of bound");
          }
          return v[j*nx+i];
        }

        const T& operator () (const int i, const int j) const {
          if((i<0) || (i>=nx)){
            ArrayExit("Array2 1st index out of bound");
          }
          if((j<0) || (j>=ny)){
            ArrayExit("Array2 2nd index out of bound");
          }
          return v[j*nx+i];
        }

        // Take data
        T* data() {return v;}
        const T* data() const {return v;}

        // Get size
        int size() const {return _size;}
    };

  // 3D Array
  template<class T>
    class Array3{

      private:
        T  *v;                // data
        int _size;            // total size
        int status;           // allocation status
        int nx = 0;           // 1st dimension
        int ny = 0;           // 2nd dimension
        int nz = 0;           // 3rd dimension

      public:

        // Constructor and deconstructor
        Array3(): v(NULL), _size(0), status(unallocated) {}

        Array3(const int nx0, const int ny0, const int nz0): 
          v(NULL), _size(0), status(unallocated){
            allocate(nx0, ny0, nz0);
          }

        Array3(const Array3<T> &b): v(NULL), _size(0), status(unallocated){
          allocate(b.nx, b.ny, b.nz);
          for(int i=0; i<b._size; i++){
            v[i] = b.v[i];
          }
        }

        ~Array3(){
          if(allocated()){ 
            deallocate();
          }
        }

        // Allocation and deallocation
        void allocate(const int nx0, const int ny0, const int nz0){
          if(allocated()){
            ArrayExit("Try to allocate an allocated Array3");
          }
          if((nx0<=0) || (ny0<=0) || (nz0<=0)){
            ArrayExit("Try to allocate non-positive size in Array3");
          }
          nx = nx0;
          ny = ny0;
          nz = nz0;
          _size = nx * ny * nz;
          v = new T[_size];
          for(int i=0; i<_size; i++){
            v[i] = 0;
          }
          status = 1;
        }

        void deallocate(){
          if(!allocated()){
            ArrayExit("Try to deallocate a unallocated Array 3");
          }
          delete [] v;
          _size = 0;
          nx = 0;
          ny = 0;
          nz = 0;
          status = unallocated;
        }

        bool allocated(){return status && (_size != 0);}  

        // Enquiry
        int Nx() const {return nx;}
        int Ny() const {return ny;}
        int Nz() const {return nz;}

        // Operator overwriting
        Array3<T>& operator = (const Array3<T>& b){
          if(allocated()) {
            if((nx != b.nx) || (ny != b.ny) || (nz != b.nz)){
              deallocate();
              allocate(b.nx, b.ny, b.nz);
            }
          }
          else {
            allocate(b.nx, b.ny, b.nz);
          }
          for(int i=0; i<b._size; i++){
            v[i] = b.v[i];
          }
          return *this;
        }

        Array3<T>& operator = (const T b){
          if(!allocated()){
            ArrayExit("Try to assign a value to a unallocated Array3");
          }
          for(int i=0; i<_size; i++){
            v[i] = b;
          }
          return *this;
        }

        T& operator () (const int i, const int j, const int k) {
          if((i<0) || (i>=nx)){
            ArrayExit("Array3 1st index out of bound");
          }
          if((j<0) || (j>=ny)){
            ArrayExit("Array3 2nd index out of bound");
          }
          if((k<0) || (k>=nz)){
            ArrayExit("Array3 3rd index out of bound");
          }
          return v[k*nx*ny+j*nx+i];
        }

        const T& operator () (const int i, const int j, const int k) const {
          if((i<0) || (i>=nx)){ 
            ArrayExit("Array3 1st index out of bound");
          }
          if((j<0) || (j>=ny)){
            ArrayExit("Array3 2nd index out of bound");
          }
          if((k<0) || (k>=nz)){
            ArrayExit("Array3 3rd index out of bound");
          }
          return v[k*nx*ny+j*nx+i];
        }

        // Take data
        T* data() {return v;}
        const T* data() const {return v;}

        // Get size
        int size() const {return _size;}
    };

};

// User defined types
typedef Array::Array1<int>    Array1i;
typedef Array::Array1<float>  Array1f;
typedef Array::Array1<double> Array1d;
typedef Array::Array1<bool>   Array1b;

typedef Array::Array2<int>    Array2i;
typedef Array::Array2<float>  Array2f;
typedef Array::Array2<double> Array2d;
typedef Array::Array2<bool>   Array2b;

typedef Array::Array3<int>    Array3i;
typedef Array::Array3<float>  Array3f;
typedef Array::Array3<double> Array3d;
typedef Array::Array3<bool>   Array3b;

#endif // ARRAY_H
