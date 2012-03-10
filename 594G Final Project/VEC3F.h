#ifndef VEC3F_H
#define VEC3F_H

//////////////////////////////////////////////////////////////////////
// This code is based on code from the (excellent) libgfx library by
// Michael Garland:
//
// http://mgarland.org/software/libgfx.html
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

class VEC3F {

public:
  // Standard constructors
  VEC3F(float x = 0, float y = 0, float z = 0) { _element[0] = x; _element[1] = y; _element[2] = z; }

  // Copy constructors & assignment operators
  VEC3F(const VEC3F& v) { *this = v; }
  VEC3F& operator=(float s) { _element[0] = _element[1] = _element[2] = s; return *this; }

  // Access methods
  operator const float*() const { return _element; }
  float& operator[](int i) { return _element[i]; }

  VEC3F& operator+=(const VEC3F& v)     { x += v.x; y += v.y; z += v.z; return *this; }
  VEC3F& operator-=(const VEC3F& v)     { x -= v.x; y -= v.y; z -= v.z; return *this; }
  VEC3F& operator/=(const VEC3F& v)     { x /= v.x; y /= v.y; z /= v.z; return *this; }
  VEC3F& operator/=(const float b)     { x /= b; y /= b; z /= b; return *this; }
  VEC3F& operator*=(const float b)     { x *= b; y *= b; z *= b;   return *this; }
  VEC3F operator+(const VEC3F &b) const { return VEC3F(x+b.x,y+b.y,z+b.z); }
  VEC3F operator-(const VEC3F &b) const { return VEC3F(x-b.x,y-b.y,z-b.z); }
  VEC3F operator*(const float b) const { return VEC3F(x*b,y*b,z*b); }
  VEC3F operator/(const float b) const { return VEC3F(x/b,y/b,z/b); }
  
  // this computes the cross product
  VEC3F operator^(const VEC3F &v) const { return VEC3F(y*v.z-v.y*z,-x*v.z+v.x*z,x*v.y-v.x*y); }
  
  // these are *element-by-element* multiplies, not dot products
  VEC3F operator*(const VEC3F& v) const { return VEC3F(x*v.x,y*v.y,z*v.z); }
  VEC3F& operator*=(const VEC3F& v)     { x *= v.x; y *= v.y; z *= v.z; return *this; }

  // *this one* does the dot product
  float dot(const VEC3F &b) const      { return x*b.x+y*b.y+z*b.z; };
  void clear() { float zero = 0.0; _element[0] = zero; _element[1] = zero; _element[2] = zero; };

  float magnitude() { return sqrt(_element[0] * _element[0] + _element[1] * _element[1] + _element[2] * _element[2]); };
  VEC3F normalize() {
    float l = _element[0] * _element[0] + 
              _element[1] * _element[1] + 
              _element[2] * _element[2];
    if (l != 1.0 && l != 0.0)
    {
      float inv = 1.0 / sqrt(l);
     _element[0] *= inv; 
     _element[1] *= inv; 
     _element[2] *= inv; 
    }
    return *this;
  };
  VEC3F normal() {
    VEC3F a = *this;
    a.normalize();
    return a;
  }

  float maxVal() { return x > y && x > z ? x : y > z ? y : z; }
  
  // the data
  union {
     struct { float x,y,z; };
     struct { float r,g,b; };
     float _element[3];
  };
};

inline std::istream &operator>>(std::istream &in, VEC3F& v)
{ return in >> v[0] >> v[1] >> v[2]; }

inline std::ostream &operator<<(std::ostream &out, VEC3F& v)
{ return out << v[0] << " " << v[1] << " " << v[2]; }

inline VEC3F operator*(const float a, const VEC3F& b)
{ return VEC3F(b.x * a, b.y * a, b.z * a); }

#endif
