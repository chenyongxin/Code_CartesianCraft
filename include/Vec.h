// Vec.h 
// Class template for vector with 3 entries.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef VEC_H
#define VEC_H

#include <cmath>
#include <iostream>

// Class template for vector
template<class T>
class Vec3{

  protected:
    T v[3];          // vector with 3 entries

  public:

    // Constructors
    Vec3(){
      for(int i=0; i<3; i++){
        v[i] = 0;
      }
    }

    Vec3(const T a, const T b, const T c){
      v[0] = a;
      v[1] = b;
      v[2] = c;
    }

    Vec3(const Vec3<T> &b){
      for(int i=0; i<3; i++){
        v[i] = b.v[i];
      }
    }

    Vec3(const T b){
      for(int i=0; i<3; i++){
        v[i] = b;
      }
    }

    // Calculation
    // Dot product
    T dot(const Vec3<T>& b) const {
      T a = 0;
      for(int i=0; i<3; i++){
        a += v[i]*b.v[i];
      }
      return a;
    }

    // Product of 3 entries
    T product() const { return v[0] * v[1] * v[2]; }

    // Cross product
    Vec3<T> cross(const Vec3<T>& b) const {
      return Vec3<T>(v[1]*b.v[2] - v[2]*b.v[1],
                     v[2]*b.v[0] - v[0]*b.v[2],
                     v[0]*b.v[1] - v[1]*b.v[0]);
    }

    // Normalise itself
    Vec3<T>& normal(){
      T a = norm2();
      for(int i=0; i<3; i++){
        v[i] /= a;
      }
      return *this;
    }

    // Get 2nd norm
    T norm2() const {
      return sqrt(dot(*this));
    }                           

    // Operator overloading
    // All entries equal to a constant
    Vec3<T>& operator = (const T b){
      for(int i=0; i<3; i++){
        v[i] = b;
      }
      return *this;
    }

    // Copy from another vector
    Vec3<T>& operator = (const Vec3<T> &b){
      for(int i=0; i<3; i++){
        v[i] = b.v[i];
      }
      return *this;
    }

    // Equivalent to another vector
    bool operator == (const Vec3<T> &b){
      bool a = true;
      for(int i=0; i<3; i++) {
        a = a && (v[i] == b.v[i]);
      }
      return a;
    }

    // Get ith reference
    T& operator () (const int i){
      if(i<0 || i>=3){
        std::cerr << std::endl << "Index out of bound in vector!" << std::endl;
        exit(1);
      }
      return v[i];
    }

    // Get ith value
    const T& operator () (const int i) const {
      if(i<0 || i>=3){
        std::cerr << std::endl << "Index out of bound in vector!" << std::endl;
        exit(1);
      }
      return v[i];
    }

    Vec3<T> operator + (const T b) const {
      return Vec3<T>(v[0]+b, v[1]+b, v[2]+b);
    }

    Vec3<T> operator + (const Vec3<T>& b) const {
      return Vec3<T>(v[0]+b.v[0], v[1]+b.v[1], v[2]+b.v[2]);
    }

    Vec3<T>& operator += (const T b){
      for(int i=0; i<3; i++){
        v[i] += b;
      }
      return *this;
    }

    Vec3<T>& operator += (const Vec3<T> &b){
      for(int i=0; i<3; i++){
        v[i] += b.v[i];
      }
      return *this;
    }

    Vec3<T> operator - (const T b) const {
      return Vec3<T>(v[0]-b, v[1]-b, v[2]-b);
    }

    Vec3<T> operator - (const Vec3<T> &b) const {
      return Vec3<T>(v[0]-b.v[0], v[1]-b.v[1], v[2]-b.v[2]);
    }

    Vec3<T>& operator -= (const T b){
      for(int i=0; i<3; i++){
        v[i] -= b;
      }
      return *this;
    }

    Vec3<T>& operator -= (const Vec3<T> &b){
      for(int i=0; i<3; i++){
        v[i] -= b.v[i];
      }
      return *this;
    }

    Vec3<T> operator * (const T b) const {
      return Vec3<T>(v[0]*b, v[1]*b, v[2]*b);
    }

    Vec3<T> operator * (const Vec3<T> &b) const {
      return Vec3<T>(v[0]*b.v[0], v[1]*b.v[1], v[2]*b.v[2]);
    }

    Vec3<T>& operator *= (const T b){
      for(int i=0; i<3; i++){
        v[i] *= b;
      }
      return *this;
    }

    Vec3<T>& operator *= (const Vec3<T> &b){
      for(int i=0; i<3; i++){
        v[i] *= b.v[i];
      }
      return *this;
    }

    Vec3<T> operator / (const T b) const {
      return Vec3<T>(v[0]/b, v[1]/b, v[2]/b);
    }

    Vec3<T> operator / (const Vec3<T> &b) const {
      return Vec3<T>(v[0]/b.v[0], v[1]/b.v[1], v[2]/b.v[2]);
    }

    Vec3<T>& operator /= (const T b){
      for(int i=0; i<3; i++){
        v[i] /= b;
      }
      return *this;
    }

    Vec3<T>& operator /= (const Vec3<T> &b){
      for(int i=0; i<3; i++){
        v[i] /= b.v[i];
      }
      return *this;
    }

    void print(){
      for(int i=0; i<3; i++){
        std::cout << v[i] << " ";
      }
      std::cout << std::endl;
    }

};

// User defined types
typedef Vec3<int>    Vec3i;
typedef Vec3<float>  Vec3f;
typedef Vec3<double> Vec3d;
typedef Vec3<bool>   Vec3b;

#endif // VEC_H
