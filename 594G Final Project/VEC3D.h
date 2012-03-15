#ifndef VEC3D_H
#define VEC3D_H

//////////////////////////////////////////////////////////////////////
// This code is based on code from the (excellent) libgfx library by
// Michael Garland:
//
// http://mgarland.org/software/libgfx.html
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

class VEC3D {
  
public:
  // Standard constructors
  VEC3D(double x = 0, double y = 0, double z = 0) { _element[0] = x; _element[1] = y; _element[2] = z; }
  
  // Copy constructors & assignment operators
  VEC3D(const VEC3D& v) { *this = v; }
  VEC3D& operator=(double s) { _element[0] = _element[1] = _element[2] = s; return *this; }
  
  // Access methods
  operator const double*() const { return _element; }
  double& operator[](int i) { return _element[i]; }
  
  VEC3D& operator+=(const VEC3D& v)     { x += v.x; y += v.y; z += v.z; return *this; }
  VEC3D& operator-=(const VEC3D& v)     { x -= v.x; y -= v.y; z -= v.z; return *this; }
  VEC3D& operator/=(const VEC3D& v)     { x /= v.x; y /= v.y; z /= v.z; return *this; }
  VEC3D& operator/=(const double b)     { x /= b; y /= b; z /= b; return *this; }
  VEC3D& operator*=(const double b)     { x *= b; y *= b; z *= b;   return *this; }
  VEC3D operator+(const VEC3D &b) const { return VEC3D(x+b.x,y+b.y,z+b.z); }
  VEC3D operator-(const VEC3D &b) const { return VEC3D(x-b.x,y-b.y,z-b.z); }
  VEC3D operator*(const double b) const { return VEC3D(x*b,y*b,z*b); }
  VEC3D operator/(const double b) const { return VEC3D(x/b,y/b,z/b); }
  
  // this computes the cross product
  VEC3D operator^(const VEC3D &v) const { return VEC3D(y*v.z-v.y*z,-x*v.z+v.x*z,x*v.y-v.x*y); }
  
  // these are *element-by-element* multiplies, not dot products
  VEC3D operator*(const VEC3D& v) const { return VEC3D(x*v.x,y*v.y,z*v.z); }
  VEC3D& operator*=(const VEC3D& v)     { x *= v.x; y *= v.y; z *= v.z; return *this; }
  
  // *this one* does the dot product
  double dot(const VEC3D &b) const      { return x*b.x+y*b.y+z*b.z; };
  void clear() { double zero = 0.0; _element[0] = zero; _element[1] = zero; _element[2] = zero; };
  
  double magnitude() { return sqrt(_element[0] * _element[0] + _element[1] * _element[1] + _element[2] * _element[2]); };
  VEC3D normalize() {
    double l = _element[0] * _element[0] + 
    _element[1] * _element[1] + 
    _element[2] * _element[2];
    if (l != 1.0 && l != 0.0)
    {
      double inv = 1.0 / sqrt(l);
      _element[0] *= inv; 
      _element[1] *= inv; 
      _element[2] *= inv; 
    }
    return *this;
  };
  VEC3D normal() {
    VEC3D a = *this;
    a.normalize();
    return a;
  }
  
  double maxVal() { return x > y && x > z ? x : y > z ? y : z; }
  
  // the data
  union {
    struct { double x,y,z; };
    struct { double r,g,b; };
    double _element[3];
  };
};

inline std::istream &operator>>(std::istream &in, VEC3D& v)
{ return in >> v[0] >> v[1] >> v[2]; }

inline std::ostream &operator<<(std::ostream &out, VEC3D& v)
{ return out << v[0] << " " << v[1] << " " << v[2]; }

inline VEC3D operator*(const double a, const VEC3D& b)
{ return VEC3D(b.x * a, b.y * a, b.z * a); }

#endif
