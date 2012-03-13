#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H


#include "PARTICLE.h"
#include "WALL.h"
#include "SPRING.h"
#include <vector>
//#include <tr1/tuple>
//#include <map>
#include "FIELD_3D.h"

#include "HW1.h"

#define h 0.045

#define PARTICLE_K 20.0 // 461.5  // Nm/kg is gas constant of water vapor
#define PARTICLE_REST_DENSITY 1000.0 // kg/m^3 is rest density of water particle
#define PARTICLE_MASS 0.012 // kg
#define PARTICLE_MEW 5.0 // 0.00089 // Ns/m^2 or Pa*s viscosity of water
#define N_SURFACE_THRESHOLD 20.0
#define PARTICLE_THETA 0.073 // N/m  coefficient of surface tension for water

#define WALL_K 10000.0 // wall spring constant
#define WALL_DAMPING -0.9 // wall damping constant

#define BOX_SIZE 0.3

using namespace std;

class PARTICLE_SYSTEM {
  
public:
  PARTICLE_SYSTEM();
  ~PARTICLE_SYSTEM();

  void updateGrid();
  
  // draw to OGL
  void draw();
  
  void stepVerlet(float dt);
  
  void stepVerletBrute(float dt);
    
  void calculateAcceleration();
  
  void calculateAccelerationBrute();
  
  void collisionForce(PARTICLE& particle, VEC3F& f_collision);
  
  void getNeighborParticles(vector<PARTICLE>& totalNeighborParticles, int x, int y, int z);
  
  float Wpoly6(float radiusSquared);
  
  void Wpoly6Gradient(VEC3F& diffPosition, float radiusSquared, VEC3F& gradient);
  
  float Wpoly6Laplacian(float radiusSquared); 
  
  float WspikyGradient(VEC3F& r); // gradient of spiky kernel
  
  float WviscosityLaplacian(float radiusSquared);
  
  //typedef std::tr1::tuple<int,int,int> gridKey;  
  //std::map<gridKey, std::vector<PARTICLE> > grid;
  
  
  FIELD_3D* grid;

private:
  // list of particles, walls, and springs being simulated
  vector<PARTICLE> _particles;
  vector<WALL>     _walls;
  vector<SPRING>   _springs;

  unsigned int _particleCount;

  // add a spring between the two particles
  void addSpring(int left, int right);
};

#endif
