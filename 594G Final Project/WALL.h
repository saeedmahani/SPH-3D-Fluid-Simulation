#ifndef WALL_H
#define WALL_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "PARTICLE.h"
#include "VEC3F.h"

using namespace std;

class WALL {
public:
  WALL(const VEC3F& normal, const VEC3F& point);

  // draw to OGL
  void draw();

  // accessors
  VEC3F& normal() { return _normal; };
  VEC3F& point()  { return _point; };

private:  
  VEC3F _normal;
  VEC3F _point;
};

#endif
