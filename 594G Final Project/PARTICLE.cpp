#include "PARTICLE.h"

VEC3F red(1,0,0);
VEC3F blue(0,0,1); 
VEC3F black(0,0,0);

int pCount = 0;

#define PARTICLE_DRAW_RADIUS 0.01//0.01 //0.006

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
PARTICLE::PARTICLE() 
{
  //_velocity[0] = 20.0; 
  //std::cout << "particle constructor\n";
}

PARTICLE::PARTICLE(const VEC3F& position) :
  _position(position)
{
  _id = pCount++;
  //_velocity[0] = 1.0;
  //std::cout << "particle constructor with vars\n";
}

PARTICLE::PARTICLE(const VEC3F& position, const VEC3F& velocity) :
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
  
  
  if (_flag)
    glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
  else
    glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
  
  glPushMatrix();
    glTranslatef(_position[0], _position[1], _position[2]);
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

