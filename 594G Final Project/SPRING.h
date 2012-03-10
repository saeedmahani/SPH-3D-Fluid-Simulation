#ifndef SPRING_H
#define SPRING_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "VEC3F.h"

class SPRING {
public:
  SPRING(int leftIndex, int rightIndex, float restLength) :
    _leftIndex(leftIndex), _rightIndex(rightIndex), _restLength(restLength) { };

  float& restLength() { return _restLength; };
  int& leftIndex()    { return _leftIndex; };
  int& rightIndex()   { return _rightIndex; };

private:
  // vertex indices the spring is attached to
  int _leftIndex;
  int _rightIndex;

  // rest length of the spring
  float _restLength;
};

#endif
