#include "PARTICLE.h"

VEC3F red(1,0,0);
VEC3F blue(0,0,1); 
VEC3F black(0,0,0);
VEC3F green(0,1,0);
VEC3F lightBlueColor(0.01,0.25,1.0);
VEC3F purpleColor(0.88,0.08,0.88);

int count = 0;

#define PARTICLE_DRAW_RADIUS 0.015//0.01 //0.006

bool PARTICLE::isSurfaceVisible = false;
bool PARTICLE::showArrows = false;
unsigned int PARTICLE::count = 0;

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////

PARTICLE::PARTICLE() 
{
}

PARTICLE::PARTICLE(const VEC3D& position) :
  _position(position)
{
  myQuadric = NULL;
  _id = count++;
}

PARTICLE::PARTICLE(const VEC3D& position, const VEC3D& velocity) :
_position(position), _velocity(velocity)
{
  myQuadric = NULL;
  _id = count++;
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
  
    
    if (showArrows) {
      // scale
      
      //glColor3f(0.2f,0.3f,0.6f);
      if (!myQuadric) {
        myQuadric = gluNewQuadric();
        gluQuadricDrawStyle(myQuadric, GLU_FILL); 
        gluQuadricNormals(myQuadric, GLU_SMOOTH);
      }
      
      
      
      double angle1 = asin(_velocity[0]) * 180.0 / M_PI;
      double angle2 = asin(_velocity[1]) * 180.0 / M_PI;
      //double angle3 = asin(_velocity[2]) * 180.0 / M_PI;
      
      
      glRotatef(-angle1, 0, 1, 0);
      glRotatef(-angle2, 1, 0, 0);
      //glRotatef(-angle3, 0, 0, 1);
      
      gluCylinder(myQuadric, 0.001, 0.001, 0.01, 10, 10);
      glTranslated(0.00, 0.01, 0.00);
      glutSolidCone(0.003, 0.01, 10, 10);
      
      glFlush();
      
      
      
      //;
    }
    else {
      glutSolidSphere(PARTICLE_DRAW_RADIUS, 10, 10);
    }
  
  glPopMatrix();
}

void PARTICLE::clearParameters() {
  _position = VEC3D();
  _velocity = VEC3D();
  _acceleration = VEC3D();
  _density = 0.0;
  _pressure = 0.0;
  
  
}

