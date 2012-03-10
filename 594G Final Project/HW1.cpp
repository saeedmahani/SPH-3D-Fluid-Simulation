///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of Example 3-1 
// from the OpenGL Red Book:
// 
// http://www.glprogramming.com/red/chapter03.html
//
///////////////////////////////////////////////////////////////////////////////


#include <cstdlib>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <glvu.h>

#include "WALL.h"
#include "PARTICLE.h"
#include "PARTICLE_SYSTEM.h"
#include <vector>

// GUI interaction stuff
GLVU glvu;

PARTICLE_SYSTEM particleSystem;

float dt = 1.0 / 4000.0;
bool animate = false;

///////////////////////////////////////////////////////////////////////
// draw coordinate axes
///////////////////////////////////////////////////////////////////////
void drawAxes()
{
  //glDisable(GL_COLOR_MATERIAL);
  // draw coordinate axes
  glPushMatrix();
  glTranslatef(-0.1f, -0.1f, -0.1f);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
  // x axis is red
  glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(10.0f, 0.0f, 0.0f, 0.0f);
  glVertex3f(1.0f, 0.0f, 0.0f);
  
  // y axis is green 
  glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(0.0f, 10.0f, 0.0f, 0.0f);
  glVertex3f(0.0f, 1.0f, 0.0f);
  
  // z axis is blue
  glColor4f(0.0f, 0.0f, 10.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(0.0f, 0.0f, 10.0f, 0.0f);
  glVertex3f(0.0f, 0.0f, 1.0f);
  glEnd();
  glLineWidth(1.0f);
  glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////
// The drawing function
///////////////////////////////////////////////////////////////////////////////
void displayCallback()
{
  glvu.BeginFrame();
  
  // clear away the previous frame
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

    
  //drawAxes();
  
  // draw the particle system and walls
  particleSystem.draw();

  // swap the buffers
  glutSwapBuffers();
  
  glvu.EndFrame();
}

///////////////////////////////////////////////////////////////////////////////
// The projection function
///////////////////////////////////////////////////////////////////////////////
void reshapeCallback(int width, int height)
{
  // set the viewport resolution (w x h)
  glViewport(0, 0, (GLsizei) width, (GLsizei) height);

  // Make ensuing transforms affect the projection matrix
  glMatrixMode(GL_PROJECTION);

  // set the projection matrix to an orthographic view
  glLoadIdentity();
  gluPerspective(65.0, (float)width / height, 0.01, 1000.0);

  // set the matric mode back to modelview
  glMatrixMode(GL_MODELVIEW);

  // set the lookat transform
  glLoadIdentity();
  gluLookAt(0, 0, 0.5, 0, 0, 0, 0, 1, 0);
}

///////////////////////////////////////////////////////////////////////////////
// Keyboard command processing function
///////////////////////////////////////////////////////////////////////////////
void keyboardCallback(unsigned char key, int x, int y)
{
  switch (key)
  {
    // quit entirely
    case 'q':
    case 'Q':
      exit(0);
      break;

    case 'a':
      animate = !animate;
      break;

  }
  glvu.Keyboard(key, x, y);
  
  glutPostRedisplay();
  
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y)
{
  glvu.Mouse(button,state,x,y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y)
{
  glvu.Motion(x,y);
}

///////////////////////////////////////////////////////////////////////////////
// Idle command processing function
///////////////////////////////////////////////////////////////////////////////
void idleCallback()
{
  if (!animate) return;

  //particleSystem.stepEuler(dt);
  //particleSystem.stepMidpoint(dt);
  //particleSystem.stepRK4(dt);
  particleSystem.stepVerlet(dt);

  glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  char title[] = "sph";
  // initialize GLUT
  glutInit(&argc, argv);
//  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE| GLUT_RGBA);
//  glutInitWindowSize(500, 500); 
//  glutInitWindowPosition(100, 100);
//  glutCreateWindow(argv[0]);
  glvu.Init(title,
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            0, 0, 800, 800);
  glShadeModel(GL_SMOOTH);

  // point GLUT to our callback functions
  glutDisplayFunc(displayCallback); 
  glutIdleFunc(idleCallback); 
  //glutReshapeFunc(reshapeCallback);
  glutKeyboardFunc(keyboardCallback);
  glutMouseFunc(glutMouseClick);
  glutMotionFunc(glutMouseMotion);
  

  // set background to black
  glClearColor(1.0, 1.0, 1.0, 1.0);

  // enable lights
  GLfloat ambient[] = {0.5,0.5,0.5};
  GLfloat diffuse[] = {0.5, 0.5, 0.5};
  GLfloat specular[] = {0.5, 0.5, 0.5};

  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  
  glvuVec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
  Eye(0, 0, 0.5),  LookAtCntr(0, 0, 0),  Up(0, 1, 0);
  
  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye, LookAtCntr, Up, Yfov, Aspect, Near, Far);
  
  glvuVec3f center(0.0, 0.0, 0.0);
  glvu.SetWorldCenter(center);
  

  //particleSystem.addWall(WALL(VEC3F(1,0,0), VEC3F(-10,0,0)));
  //particleSystem.addWall(WALL(VEC3F(-1,0,0), VEC3F(10,0,0)));
  //particleSystem.addWall(WALL(VEC3F(0,1,0), VEC3F(0,-10,0)));

  // Let GLUT take over
  glutMainLoop();
  
  return 0;
}
