#include "PARTICLE_SYSTEM.h"


///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
PARTICLE_SYSTEM::PARTICLE_SYSTEM() : _particleCount(0)
{
  int gridSize = (int)ceil(BOX_SIZE/h);
  grid = new FIELD_3D(gridSize, gridSize, gridSize);

  vector<PARTICLE>& firstGridCell = (*grid)(0,0,0);
  firstGridCell.reserve(200);
  
  
#ifdef BRUTE
  cout << "using BRUTE neighbor method" << endl;
#else
  cout << "using GRID neighbor method" << endl;
#endif
  

  bool offset = false;
  
  for (float y = 0; y < BOX_SIZE/2.0 - 0.01; y+= h/2.0) {
    for (float x = BOX_SIZE/4.0; x < BOX_SIZE/2.0 - 0.01; x += h/2.0) {
      for (float z = -BOX_SIZE/2.0; z < BOX_SIZE/2.0 - 0.01; z+= h/2.0) {
#ifdef BRUTE
         _particles.push_back(PARTICLE(VEC3F((offset ? -h/2.0 : 0) + x, y,z)));
#else
        firstGridCell.push_back(PARTICLE(VEC3F((offset ? -h/2.0 : 0) + x, y,z)));
#endif
        _particleCount++;
        //cout << firstGridCell.size() << endl;
         
      }
    }
    offset = !offset;
  }
  
#ifndef BRUTE
  updateGrid(); // to properly position all the particles initially
#endif
   
  printf("simulating %d particles\n", _particleCount);
  
  _walls.push_back(WALL(VEC3F(0,0,1), VEC3F(0,0,-BOX_SIZE/2.0)));
  _walls.push_back(WALL(VEC3F(1,0,0), VEC3F(-BOX_SIZE/2.0,0,0)));
  _walls.push_back(WALL(VEC3F(-1,0,0), VEC3F(BOX_SIZE/2.0,0,0)));
  _walls.push_back(WALL(VEC3F(0,1,0), VEC3F(0,-BOX_SIZE/2.0,0)));
  _walls.push_back(WALL(VEC3F(0,0,-1), VEC3F(0,0,BOX_SIZE/2.0)));
  //_walls.push_back(WALL(VEC3F(0,-1,0), VEC3F(0,BOX_SIZE/2.0,0)));  // top wall

  
}


PARTICLE_SYSTEM::~PARTICLE_SYSTEM()
{
  if (grid) delete grid;
}

// to update the grid cells particles are located in
// should be called right after particle positions are updated
void PARTICLE_SYSTEM::updateGrid() {
  
  //static float invGridSpacing = 1.0/h;
  
  for (unsigned int x = 0; x < (*grid).xRes(); x++) {
    for (unsigned int y = 0; y < (*grid).yRes(); y++) {
      for (unsigned int z = 0; z < (*grid).zRes(); z++) {
        
        vector<PARTICLE>& particles = (*grid)(x,y,z);
        
        //cout << particles.size() << "p's in this grid" << endl;
                
        for (int p = 0; p < particles.size(); p++) {
          
          PARTICLE& particle = particles[p];
          
          int newGridCellX = (int)floor((particle.position().x+BOX_SIZE/2.0)/h); 
          int newGridCellY = (int)floor((particle.position().y+BOX_SIZE/2.0)/h);
          int newGridCellZ = (int)floor((particle.position().z+BOX_SIZE/2.0)/h);
          
          //cout << "particle position: " << particle.position() << endl;
          //cout << "particle cell pos: " << newGridCellX << " " << newGridCellY << " " << newGridCellZ << endl;
        
          if (newGridCellX < 0)
            newGridCellX = 0;
          else if (newGridCellX >= (*grid).xRes())
            newGridCellX = (*grid).xRes() - 1;
          if (newGridCellY < 0)
            newGridCellY = 0;
          else if (newGridCellY >= (*grid).yRes())
            newGridCellY = (*grid).yRes() - 1;
          if (newGridCellZ < 0)
            newGridCellZ = 0;
          else if (newGridCellZ >= (*grid).zRes())
            newGridCellZ = (*grid).zRes() - 1;
          
          //cout << "particle cell pos: " << newGridCellX << " " << newGridCellY << " " << newGridCellZ << endl;

          
          // check if particle has moved
          
          if (x != newGridCellX || y != newGridCellY || z != newGridCellZ) {
            
            // move the particle to the new grid cell
            
            (*grid)(newGridCellX, newGridCellY, newGridCellZ).push_back(particle);
            
            // remove it from it's previous grid cell
            
            particles[p] = particles.back();
            particles.pop_back();
            p--; // important! make sure to redo this index, since a new particle will (probably) be there
          }
          
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// OGL drawing
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::draw() 
{ 
  VEC3F black(0,0,0); 
  VEC3F blue(0,0,1); 
  VEC3F white(1,1,1); 
  VEC3F green(34.0 / 255, 139.0 / 255, 34.0 / 255);
  float shininess = 20.0;

  // draw the particles
  glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
  glMaterialfv(GL_FRONT, GL_SPECULAR, white);
  glMaterialfv(GL_FRONT, GL_SHININESS, &shininess);
  
  //for (unsigned int x = 0; x < _particles.size(); x++)
  //  _particles[x].draw();
  
#ifdef BRUTE
  
  for (unsigned int x = 0; x < _particles.size(); x++)
    _particles[x].draw();
    
#else
  
  for (int gridCellIndex = 0; gridCellIndex < (*grid).cellCount(); gridCellIndex++) {
    
    vector<PARTICLE>& particles = (*grid).data()[gridCellIndex];
    
    for (int p = 0; p < particles.size(); p++) {
      
      PARTICLE& particle = particles[p];
      
      particle.draw();
      
    }
  }

#endif

  
  glutWireCube(BOX_SIZE);
  
//  glPushMatrix();
//  glTranslatef(BOX_SIZE/2.0, 0, 0);
//  glutWireSphere(0.01, 10, 10);
//  glPopMatrix();
  
  /*
  
  // draw the springs
  glLineWidth(5.0);
  glBegin(GL_LINES);
  for (unsigned int x = 0; x < _springs.size(); x++)
  {
    int left = _springs[x].leftIndex();
    int right = _springs[x].rightIndex();
    VEC3F& p0 = _particles[left].position();
    VEC3F& p1 = _particles[right].position();
    glVertex3fv(p0);
    glVertex3fv(p1);
  }
  glEnd();

  // draw the walls
  glMaterialfv(GL_FRONT, GL_SPECULAR, green);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
  for (unsigned int x = 0; x < _walls.size(); x++)
    _walls[x].draw();
   */
   
}

///////////////////////////////////////////////////////////////////////////////
// Verlet integration
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::stepVerlet(float dt)
{
  
  calculateAcceleration();
  
  for (unsigned int gridCellIndex = 0; gridCellIndex < (*grid).cellCount(); gridCellIndex++) {
    
    vector<PARTICLE>& particles = (*grid).data()[gridCellIndex];
    
    for (unsigned int p = 0; p < particles.size(); p++) {
      
      PARTICLE& particle = particles[p];
      
      VEC3F newPosition = particle.position() + particle.velocity()*dt + particle.acceleration()*dt*dt;
      VEC3F newVelocity = (newPosition - particle.position()) / dt;
      
      particle.position() = newPosition;
      particle.velocity() = newVelocity;
    }
  }
}

float totalDensity = 0.0;
VEC3F totalPressure;
int neighborsVisited = 0;

void PARTICLE_SYSTEM::calculateAcceleration() {
    
  // to conaint a vector of neighbors belonging to each particle (same index as _particles)
  //vector<PARTICLE> **neighborArray = new vector<PARTICLE>*[_particleCount]; // make this vector instead?
  
  
  ///////////////////
  // STEP 1: UPDATE DENSITY & PRESSURE OF EACH PARTICLE
  
  //cout << "gridCellCount: " << (*grid).cellCount() << endl;
  //cout << "grid xRes: " << (*grid).xRes() << endl;
  
  for (int x = 0; x < (*grid).xRes(); x++) {
    for (int y = 0; y < (*grid).yRes(); y++) {
      for (int z = 0; z < (*grid).zRes(); z++) {
        
        vector<PARTICLE>& particles = (*grid)(x,y,z);
                
        for (int p = 0; p < particles.size(); p++) {
          
          PARTICLE& particle = particles[p];

          particle.density() = 0.0;
          
          // now iteratate through neighbors
          
          for (int offsetX = -1; offsetX <= 1; offsetX++) {
            if (x+offsetX < 0) continue;
            if (x+offsetX >= (*grid).xRes()) break;
            
            for (int offsetY = -1; offsetY <= 1; offsetY++) {
              if (y+offsetY < 0) continue;
              if (y+offsetY >= (*grid).yRes()) break;
              
              for (int offsetZ = -1; offsetZ <= 1; offsetZ++) {
                if (z+offsetZ < 0) continue;
                if (z+offsetZ >= (*grid).zRes()) break;
                
                vector<PARTICLE>& neighborGridCellParticles = (*grid)(x+offsetX, y+offsetY, z+offsetZ);
              
                for (int i = 0; i < neighborGridCellParticles.size(); i++) {
                  
                  PARTICLE& neighborParticle = neighborGridCellParticles[i];
                  
                  VEC3F diffPosition = particle.position() - neighborParticle.position();
                  
                  float radiusSquared = diffPosition.dot(diffPosition);
                  
                  if (radiusSquared <= h*h)
                    particle.density() += Wpoly6(radiusSquared);
                  
                }
              }
            }
          }
          
          particle.density() *= PARTICLE_MASS;
                        
          // p = k(density - density_rest)
          
          particle.pressure() = PARTICLE_K * (particle.density() - PARTICLE_REST_DENSITY);
          
        }
      }
    }
  }
  
  
  //////////// testing total density
  
  for (int gridCellIndex = 0; gridCellIndex < (*grid).cellCount(); gridCellIndex++) {
    
    vector<PARTICLE>& particles = (*grid).data()[gridCellIndex];
    
    for (int p = 0; p < particles.size(); p++) {
      
      PARTICLE& particle = particles[p];
      
      totalDensity += particle.density();
      
    }
  }
  
  
  //////////
  
  
  //cout << "total density (GRID): " << totalDensity << endl;
  
  ///////////////////
  // STEP 2: COMPUTE FORCES FOR ALL PARTICLES
  
  for (int x = 0; x < (*grid).xRes(); x++) {
    for (int y = 0; y < (*grid).yRes(); y++) {
      for (int z = 0; z < (*grid).zRes(); z++) {
        
        vector<PARTICLE>& particles = (*grid)(x,y,z);
        
        for (int p = 0; p < particles.size(); p++) {
          
          PARTICLE& particle = particles[p];
          
          //cout << "particle id: " << particle.id() << endl;
          
          VEC3F f_pressure, 
          f_viscosity, 
          f_surface, 
          f_gravity(0.0, particle.density() * -9.80665, 0.0),
          n, 
          colorFieldNormal,
          colorFieldLaplacian;
                    
          // now iteratate through neighbors
          
          for (int offsetX = -1; offsetX <= 1; offsetX++) {
            if (x+offsetX < 0) continue;
            if (x+offsetX >= (*grid).xRes()) break;
            
            for (int offsetY = -1; offsetY <= 1; offsetY++) {
              if (y+offsetY < 0) continue;
              if (y+offsetY >= (*grid).yRes()) break;
              
              for (int offsetZ = -1; offsetZ <= 1; offsetZ++) {
                if (z+offsetZ < 0) continue;
                if (z+offsetZ >= (*grid).zRes()) break;
                
                vector<PARTICLE>& neighborGridCellParticles = (*grid)(x+offsetX, y+offsetY, z+offsetZ);
                
                for (int i = 0; i < neighborGridCellParticles.size(); i++) {
                  
                  PARTICLE& neighbor = neighborGridCellParticles[i];
                  
                  VEC3F diffPosition = particle.position() - neighbor.position();
                  //VEC3F diffPositionNormalized = diffPosition.normal(); // need?
                  float radiusSquared = diffPosition.dot(diffPosition);
                  
                  if (radiusSquared <= h*h) {
                    
                    
                    
                    if (radiusSquared > 0.0) {
                      
                      neighborsVisited++;
                      //cout << neighborsVisited << endl;
                      //cout << neighbor.id() << "," << endl;
                      
                      VEC3F gradient;
                                            
                      Wpoly6Gradient(diffPosition, radiusSquared, gradient);
                                            
                      f_pressure += (particle.pressure() + neighbor.pressure()) / (2.0 * neighbor.density()) * gradient;
                      
                      colorFieldNormal += gradient / neighbor.density();
                    }
                    
                    f_viscosity += (neighbor.velocity() - particle.velocity()) * WviscosityLaplacian(radiusSquared) / neighbor.density();
                    
                    colorFieldLaplacian += Wpoly6Laplacian(radiusSquared) / neighbor.density();
                  }
                  
                }
              }
            }
          } // end of neighbor grid cell iteration
          
          f_pressure *= -PARTICLE_MASS;
          
          totalPressure += f_pressure;
          
          //cout << "pressure: " << f_pressure << endl;
          
          f_viscosity *= PARTICLE_MEW * PARTICLE_MASS;
          
          colorFieldNormal *= PARTICLE_MASS;
          
          colorFieldLaplacian *= PARTICLE_MASS;
          
          // ADD IN SPH FORCES
          
          particle.acceleration() = (f_pressure + f_viscosity + f_surface + f_gravity) / particle.density();
          
          //cout << particle.acceleration() << endl;
          
          // EXTERNAL FORCES HERE (USER INTERACTION, SWIRL)
          
          VEC3F f_collision; 
          collisionForce(particle, f_collision);
          
          
          
          particle.acceleration() += (f_collision) / PARTICLE_MASS;
        } // end of particle iteration
      }
    }
  }
  
  //cout << "neighbors visited: " << neighborsVisited << endl;
  //cout << "total pressure: " << totalPressure << endl;
}


void PARTICLE_SYSTEM::collisionForce(PARTICLE& particle, VEC3F& f_collision) {
    
  for (unsigned int i = 0; i < _walls.size(); i++) {
    
    WALL& wall = _walls[i];
    
    float d = (wall.point() - particle.position()).dot(wall.normal()) + 0.01; // particle radius
    
    if (d > 0.0) {
      //particle.position() += d * wall.normal();
      //particle.velocity() -= particle.velocity().dot(wall.normal()) * 1.94 * wall.normal();
      
      
      
      f_collision += WALL_K * wall.normal() * d;
      f_collision += WALL_DAMPING * particle.velocity().dot(wall.normal()) * wall.normal();
    }
  }
}


float PARTICLE_SYSTEM::Wpoly6(float radiusSquared) {  // checked
    
  static float coefficient = 315.0/(64.0*M_PI*pow(h,9));
  static float hSquared = h*h;
  
  return coefficient * pow(hSquared-radiusSquared, 3);
}

void PARTICLE_SYSTEM::Wpoly6Gradient(VEC3F& diffPosition, float radiusSquared, VEC3F& gradient) {  
    
  static float coefficient = -6.0*315.0/(64.0*M_PI*pow(h,9));
  static float hSquared = h*h;
  
  float weight = coefficient * pow(hSquared-radiusSquared, 2);
  
  gradient = weight * diffPosition;
}

float PARTICLE_SYSTEM::Wpoly6Laplacian(float radiusSquared) {
    
  static float coefficient = 315.0/(64.0*M_PI*pow(h,9));
  static float hSquared = h*h;
  
  float hSquaredMinusRadiusSquared = hSquared - radiusSquared;
  
  return coefficient * (-18.0*pow(hSquaredMinusRadiusSquared, 2) + 24.0*radiusSquared*hSquaredMinusRadiusSquared);
}

/*
float PARTICLE_SYSTEM::WspikyGradient(VEC3F& r) {  // checked
  
  float radius = r.magnitude();
  
  //cout << "spiky1 r: " << r << endl;
  //cout << "spiky1 radius: " << radius << endl;
  
  if (radius < 0 || radius > h) {
    return 0.0;
  }
  else {
    //cout << "spiky1:" << (15.0/(pi*pow(h,6)) * -3.0 * pow( h-radius, 2)) << endl;
    return 15.0/(pi*pow(h,6)) * -3.0 * pow( h-radius, 2);
  }
  
}
 */

float PARTICLE_SYSTEM::WviscosityLaplacian(float radiusSquared) {  // checked
  
  static float coefficient = 45.0/(M_PI*pow(h,6));
  
  float radius = sqrt(radiusSquared);
  
  return coefficient * (h - radius);    
}


//std::cout << grid.size() << "\n";
//gridKey key = std::tr1::make_tuple(2,3,4);
//vector<PARTICLE> *hi = &grid[key];
//if (gridKey.find("

/* cal f_surface here
 
 if (n_mag > N_SURFACE_THRESHOLD) {  
 
 VEC3F colorFieldLaplacian; // laplacian of colorfield
 
 for (int j = 0; j < numNeighborParticles; j++) {
 
 PARTICLE& neighbor = neighbors[j];
 
 VEC3F diffPosition = particle.position() - neighbor.position();
 
 colorFieldLaplacian += PARTICLE_MASS / neighbor.density() * Wpoly6Laplacian(diffPosition) * diffPosition.normal();
 
 }
 
 particle.flag() = true;
 
 f_surface = -PARTICLE_THETA * colorFieldLaplacian * (n / n_mag);
 
 }
 else {
 particle.flag() = false;
 }
 */

void PARTICLE_SYSTEM::calculateAccelerationBrute() {
    
  unsigned int numParticles = _particles.size();
  
  // to conaint a vector of neighbors belonging to each particle (same index as _particles)
//  vector<PARTICLE> **neighborArray = new vector<PARTICLE>*[numParticles]; // make this vector instead?
  
  
  ///////////////////
  // STEP 1: UPDATE DENSITY & PRESSURE OF EACH PARTICLE
  
  for (int i = 0; i < _particleCount; i++) {
    
    // grab ith particle reference
    
    PARTICLE& particle = _particles[i];
    
    // now iteratate through neighbors
    
    particle.density() = 0.0;
            
    for (int j = 0; j < _particleCount; j++) {
      
      PARTICLE& neighborParticle = _particles[j];
      
      VEC3F diffPosition = particle.position() - neighborParticle.position();
      
      float radiusSquared = diffPosition.dot(diffPosition);
      
      if (radiusSquared <= h*h)
        particle.density() += Wpoly6(radiusSquared);
    }
    
    particle.density() *= PARTICLE_MASS;
    
    // p = k(density - density_rest)
    
    particle.pressure() = PARTICLE_K * (particle.density() - PARTICLE_REST_DENSITY);
  
    totalDensity += particle.density();
  }
  
  
  //cout << "total density (BRUTE): " << totalDensity << endl;

  
  
  ///////////////////
  // STEP 2: COMPUTE FORCES FOR ALL PARTICLES
  
  for (int i = 0; i < numParticles; i++) {
    
    PARTICLE& particle = _particles[i];
    
    //cout << "particle id: " << particle.id() << endl;
    
    VEC3F f_pressure, 
    f_viscosity, 
    f_surface, 
    f_gravity(0.0, particle.density() * -9.80665, 0.0),
    n, 
    colorFieldNormal,
    colorFieldLaplacian; // n is gradient of colorfield
    //float n_mag;
    
    for (int j = 0; j < _particleCount; j++) {
      PARTICLE& neighbor = _particles[j];
      
      VEC3F diffPosition = particle.position() - neighbor.position();
      VEC3F diffPositionNormalized = diffPosition.normal(); // need?
      float radiusSquared = diffPosition.dot(diffPosition);
      
      if (radiusSquared <= h*h) {
        
                
        if (radiusSquared > 0.0) {
          
          neighborsVisited++;
          //cout << neighborsVisited << endl;
          
          //cout << neighbor.id() << endl;
          
          VEC3F gradient;
                    
          Wpoly6Gradient(diffPosition, radiusSquared, gradient);
                    
          f_pressure += (particle.pressure() + neighbor.pressure()) / (2.0 * neighbor.density()) * gradient;
          
          colorFieldNormal += gradient / neighbor.density();
        }
        
        f_viscosity += (neighbor.velocity() - particle.velocity()) * WviscosityLaplacian(radiusSquared) / neighbor.density();
        
        colorFieldLaplacian += Wpoly6Laplacian(radiusSquared) / neighbor.density();
      }
      
    }
    
    f_pressure *= -PARTICLE_MASS;
    
    totalPressure += f_pressure;
    
    f_viscosity *= PARTICLE_MEW * PARTICLE_MASS;
    
    colorFieldNormal *= PARTICLE_MASS;
    
    colorFieldLaplacian *= PARTICLE_MASS;
    
    /* cal f_surface here
     
     if (n_mag > N_SURFACE_THRESHOLD) {  
     
     VEC3F colorFieldLaplacian; // laplacian of colorfield
     
     for (int j = 0; j < numNeighborParticles; j++) {
     
     PARTICLE& neighbor = neighbors[j];
     
     VEC3F diffPosition = particle.position() - neighbor.position();
     
     colorFieldLaplacian += PARTICLE_MASS / neighbor.density() * Wpoly6Laplacian(diffPosition) * diffPosition.normal();
     
     }
     
     particle.flag() = true;
     
     f_surface = -PARTICLE_THETA * colorFieldLaplacian * (n / n_mag);
     
     }
     else {
     particle.flag() = false;
     }
     */
    
    // ADD IN SPH FORCES
    
    particle.acceleration() = (f_pressure + f_viscosity + f_surface + f_gravity) / particle.density();
    
    
    // EXTERNAL FORCES HERE (USER INTERACTION, SWIRL)
    
    VEC3F f_collision;
    
    collisionForce(particle, f_collision);    
  
    
    particle.acceleration() += (f_collision) / PARTICLE_MASS;
  }
  
  //cout << "neighbors visited: " << neighborsVisited << endl;
  //cout << "total pressure: " << totalPressure << endl;
  
  //delete[] neighborArray;
}

///////////////////////////////////////////////////////////////////////////////
// Verlet integration
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::stepVerletBrute(float dt)
{
  
   unsigned int numParticles = _particles.size();
   
   calculateAccelerationBrute();
   
   for (unsigned int i = 0; i < numParticles; i++) {
   
   PARTICLE& particle = _particles[i];
   
   VEC3F newPosition = particle.position() + particle.velocity()*dt + particle.acceleration()*dt*dt;
   VEC3F newVelocity = (newPosition - particle.position()) / dt;
   
   particle.position() = newPosition;
   particle.velocity() = newVelocity;
   
   // collision detection
   
   // moved to calculateAcceleration()
   
   }
   
}
