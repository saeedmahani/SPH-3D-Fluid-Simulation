#include "PARTICLE.h"

VEC3F red(1,0,0);
VEC3F blue(0,0,1); 
VEC3F black(0,0,0);
VEC3F green(0,1,0);
VEC3F lightBlueColor(0.01,0.25,1.0);
VEC3F purpleColor(0.88,0.08,0.88);

int pCount = 0;

#define PARTICLE_DRAW_RADIUS 0.015//0.01 //0.006

bool PARTICLE::isSurfaceVisible = false;

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
PARTICLE::PARTICLE() 
{
  //_velocity[0] = 20.0; 
  //std::cout << "particle constructor\n";
}

PARTICLE::PARTICLE(const VEC3D& position) :
  _position(position)
{
  _id = pCount++;
  //_velocity[0] = 1.0;
  //std::cout << "particle constructor with vars\n";
}

PARTICLE::PARTICLE(const VEC3D& position, const VEC3D& velocity) :
_position(position), _velocity(velocity)
{
  _id = pCount++;  
}



///////////////////////////////////////////////////////////////////////////////
// OGL drawing
///////////////////////////////////////////////////////////////////////////////
void PARTICLE::draw() 
{
  /*
  if (_flag) 
    glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
  else 
    glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
  */
  
  
  if (_flag && isSurfaceVisible)
    glMaterialfv(GL_FRONT, GL_DIFFUSE, purpleColor);
  else
    glMaterialfv(GL_FRONT, GL_DIFFUSE, lightBlueColor);
  
  glPushMatrix();
    glTranslated(_position[0], _position[1], _position[2]);
    glutSolidSphere(PARTICLE_DRAW_RADIUS, 10, 10);
  
    /*
    glColor3f(1.0, 0, 0);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
    glBegin(GL_LINES);
    
    glVertex3f(0, 0, 0);
    glVertex3fv(normal*0.001);
    
    glEnd();
  
   */
  
  glPopMatrix();
  
  
  
  
  
  
}

