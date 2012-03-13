#ifndef PARTICLE_H
#define PARTICLE_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "SPRING.h"
#include "VEC3F.h"
#include <vector>

using namespace std;



class PARTICLE {
  
  
  
public:
  
  //static unsigned int count;
  PARTICLE();
  PARTICLE(const VEC3F& position);

  // draw to OGL
  void draw();

  // clear all previous accumulated forces
  void clearForce() { _force *= 0; };

  // accumulate forces
  void addForce(VEC3F newForce) { _force += newForce; };
  
  void calculateAcceleration();

  // accessors
  VEC3F& position() { return _position; };
  VEC3F& velocity() { return _velocity; };
  VEC3F& acceleration() { return _acceleration; }
  VEC3F& force()    { return _force; };
  float& density()  { return _density; };
  float& pressure() { return _pressure; };
  bool& flag() { return _flag; };
  int& id() { return _id; };
  
private:  
  VEC3F _position;
  VEC3F _velocity;
  VEC3F _force;
  VEC3F _acceleration;
  float _density;
  float _pressure;
  bool _flag;
  int _id;
};

//unsigned int PARTICLE::count = 0;

#endif
