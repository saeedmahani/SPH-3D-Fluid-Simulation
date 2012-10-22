#include "PARTICLE_SYSTEM.h"




unsigned int iteration = 0;
int scenario;

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
PARTICLE_SYSTEM::PARTICLE_SYSTEM() : 
_isGridVisible(false), surfaceThreshold(0.01), gravityVector(0.0,GRAVITY_ACCELERATION,0.0), grid(NULL)
{
  loadScenario(INITIAL_SCENARIO);
}

void PARTICLE_SYSTEM::loadScenario(int newScenario) {
  
  scenario = newScenario;
  
  // remove all particles
  
  if (grid) {
    delete grid;
    
  }
    
  _walls.clear();
  
  
  // reset params
  
  PARTICLE::count = 0;
  
  iteration = 0;
    
  
  if (scenario == SCENARIO_DAM) {
    
    // create long grid
    
    boxSize.x = BOX_SIZE*2.0;
    boxSize.y = BOX_SIZE;
    boxSize.z = BOX_SIZE/2.0;
    
    int gridXRes = (int)ceil(boxSize.x/h);
    int gridYRes = (int)ceil(boxSize.y/h);
    int gridZRes = (int)ceil(boxSize.z/h);
    
    grid = new FIELD_3D(gridXRes, gridYRes, gridZRes);
    
    
    // add walls 
    
    _walls.push_back(WALL(VEC3D(0,0,1), VEC3D(0,0,-boxSize.z/2.0))); // back
    _walls.push_back(WALL(VEC3D(0,0,-1), VEC3D(0,0,boxSize.z/2.0))); // front
    _walls.push_back(WALL(VEC3D(1,0,0), VEC3D(-boxSize.x/2.0,0,0)));     // left
    _walls.push_back(WALL(VEC3D(-1,0,0), VEC3D(boxSize.x/2.0,0,0)));     // right
    _walls.push_back(WALL(VEC3D(0,1,0), VEC3D(0,-boxSize.y/2.0,0))); // bottom
    
    vector<PARTICLE>& firstGridCell = (*grid)(0,0,0);
    
    // add particles
    
    for (double y = -boxSize.y/2.0; y < boxSize.y/2.0; y+= h/2.0) {
      for (double x = -boxSize.x/2.0; x < -boxSize.x/4.0; x += h/2.0) {
        for (double z = -boxSize.z/2.0; z < boxSize.z/2.0; z+= h/2.0) {
          firstGridCell.push_back(PARTICLE(VEC3D(x, y,z)));
        }
      }
    }
    
    cout << "Loaded dam scenario" << endl;
    cout << "Grid size is " << (*grid).xRes() << "x" << (*grid).yRes() << "x" << (*grid).zRes() << endl;
    cout << "Simulating " << PARTICLE::count << " particles" << endl;
    
  }
  else if (scenario == SCENARIO_CUBE) {
    
    // create cubed grid
    
    boxSize.x = BOX_SIZE;
    boxSize.y = BOX_SIZE;
    boxSize.z = BOX_SIZE;
    
    int gridXRes = (int)ceil(boxSize.x/h);
    int gridYRes = (int)ceil(boxSize.y/h);
    int gridZRes = (int)ceil(boxSize.z/h);
    
    grid = new FIELD_3D(gridXRes, gridYRes, gridZRes);
    
    // add walls 
    
    _walls.push_back(WALL(VEC3D(0,0,1), VEC3D(0,0,-boxSize.z/2.0))); // back
    _walls.push_back(WALL(VEC3D(0,0,-1), VEC3D(0,0,boxSize.z/2.0))); // front
    _walls.push_back(WALL(VEC3D(1,0,0), VEC3D(-boxSize.x/2.0,0,0))); // left
    _walls.push_back(WALL(VEC3D(-1,0,0), VEC3D(boxSize.x/2.0,0,0))); // right
    _walls.push_back(WALL(VEC3D(0,1,0), VEC3D(0,-boxSize.y/2.0,0))); // bottom
    
    vector<PARTICLE>& firstGridCell = (*grid)(0,0,0);
    
    // add particles
    
    for (double y = 0; y < boxSize.y; y+= h/2.0) {
      for (double x = -boxSize.x/4.0; x < boxSize.x/4.0; x += h/2.0) {
        for (double z = -boxSize.z/4.0; z < boxSize.z/4.0; z+= h/2.0) {
          firstGridCell.push_back(PARTICLE(VEC3D(x,y,z)));
        }
      }
    }
    
    cout << "Loaded cube scenario" << endl;
    cout << "Grid size is " << (*grid).xRes() << "x" << (*grid).yRes() << "x" << (*grid).zRes() << endl;
    cout << "Simulating " << PARTICLE::count << " particles" << endl;
  }
  else if (scenario == SCENARIO_FAUCET) {
    
    // create cubed grid
    
    boxSize.x = BOX_SIZE;
    boxSize.y = BOX_SIZE;
    boxSize.z = BOX_SIZE;
    
    int gridXRes = (int)ceil(boxSize.x/h);
    int gridYRes = (int)ceil(boxSize.y/h);
    int gridZRes = (int)ceil(boxSize.z/h);
    
    grid = new FIELD_3D(gridXRes, gridYRes, gridZRes);
    
    // add walls 
    
    _walls.push_back(WALL(VEC3D(0,0,1), VEC3D(0,0,-boxSize.z/2.0))); // back
    _walls.push_back(WALL(VEC3D(0,0,-1), VEC3D(0,0,boxSize.z/2.0))); // front
    _walls.push_back(WALL(VEC3D(1,0,0), VEC3D(-boxSize.x/2.0,0,0))); // left
    _walls.push_back(WALL(VEC3D(-1,0,0), VEC3D(boxSize.x/2.0,0,0))); // right
    _walls.push_back(WALL(VEC3D(0,1,0), VEC3D(0,-boxSize.y/2.0,0))); // bottom
    
    cout << "Loaded faucet scenario" << endl;
    cout << "Grid size is " << (*grid).xRes() << "x" << (*grid).yRes() << "x" << (*grid).zRes() << endl;
  }
  
  
  updateGrid();

  
  
}

void PARTICLE_SYSTEM::generateFaucetParticleSet() {
  
  VEC3D initialVelocity(-1.8,-1.8,0);
  
  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0,BOX_SIZE+h*0.6, 0), initialVelocity);
  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE, 0), initialVelocity);
  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE+h*-0.6, 0), initialVelocity);
  
  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE+h*0.3, h*0.6), initialVelocity);
  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE+h*0.3, h*-0.6), initialVelocity);
  
  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE+h*-0.3, h*0.6), initialVelocity);
  addParticle(VEC3D(BOX_SIZE/2.0-h/2.0, BOX_SIZE+h*-0.3, h*-0.6), initialVelocity);


}

void PARTICLE_SYSTEM::addParticle(const VEC3D& position, const VEC3D& velocity) {
  
#ifdef BRUTE
  _particles.push_back(PARTICLE(position, velocity));
#else
  (*grid)(0,0,0).push_back(PARTICLE(position, velocity));
#endif
  
  //_particleCount++;
}

void PARTICLE_SYSTEM::addParticle(const VEC3D& position) {
  
  addParticle(position, VEC3D());
}

PARTICLE_SYSTEM::~PARTICLE_SYSTEM()
{
  if (grid) delete grid;
}

void PARTICLE_SYSTEM::toggleGridVisble() {
  _isGridVisible = !_isGridVisible;
}

void PARTICLE_SYSTEM::toggleSurfaceVisible() {
  PARTICLE::isSurfaceVisible = !PARTICLE::isSurfaceVisible;
}

void PARTICLE_SYSTEM::toggleTumble() {
  _tumble = !_tumble;
}

void PARTICLE_SYSTEM::toggleGravity() {
  if (gravityVector.magnitude() > 0.0) {
    gravityVector = VEC3D(0,0,0);
  }
  else {
    gravityVector = VEC3D(0,GRAVITY_ACCELERATION,0);
  }
}

void PARTICLE_SYSTEM::toggleArrows() {
  PARTICLE::showArrows = !PARTICLE::showArrows;
}

void PARTICLE_SYSTEM::setGravityVectorWithViewVector(VEC3D viewVector) {
  
  if (_tumble)
    gravityVector = viewVector * GRAVITY_ACCELERATION;
  
}

// to update the grid cells particles are located in
// should be called right after particle positions are updated
void PARTICLE_SYSTEM::updateGrid() {
    
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
  static VEC3F blackColor(0,0,0); 
  static VEC3F blueColor(0,0,1); 
  static VEC3F whiteColor(1,1,1);
  static VEC3F greyColor(0.2, 0.2, 0.2);
  static VEC3F lightGreyColor(0.8,0.8,0.8);
  //static VEC3F greenColor(34.0 / 255, 139.0 / 255, 34.0 / 255);
  static float shininess = 10.0;

  // draw the particles
  glEnable(GL_LIGHTING);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, blueColor);
  glMaterialfv(GL_FRONT, GL_SPECULAR, whiteColor);
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
  
  glDisable(GL_LIGHTING);
  
  if (_isGridVisible) {
  
    // draw the grid
    
    glColor3fv(lightGreyColor);
    
    //double offset = -BOX_SIZE/2.0+h/2.0;
    
    for (int x = 0; x < grid->xRes(); x++) {
      for (int y = 0; y < grid->yRes(); y++) {
        for (int z = 0; z < grid->zRes(); z++) {
          glPushMatrix();
          
          glTranslated(x*h-boxSize.x/2.0+h/2.0, y*h-boxSize.y/2.0+h/2.0, z*h-boxSize.z/2.0+h/2.0);
          glutWireCube(h);
          
          glPopMatrix();
        }
      }
    }
    
  }
   

#endif

  glColor3fv(greyColor);
  
  glPopMatrix();
  glScaled(boxSize.x, boxSize.y, boxSize.z);
  glutWireCube(1.0);
  glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////
// Verlet integration
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::stepVerlet(double dt)
{
  
  calculateAcceleration();
  
  for (unsigned int gridCellIndex = 0; gridCellIndex < (*grid).cellCount(); gridCellIndex++) {
    
    vector<PARTICLE>& particles = (*grid).data()[gridCellIndex];
    
    for (unsigned int p = 0; p < particles.size(); p++) {
      
      PARTICLE& particle = particles[p];
      
      VEC3D newPosition = particle.position() + particle.velocity()*dt + particle.acceleration()*dt*dt;
      VEC3D newVelocity = (newPosition - particle.position()) / dt;
      
      particle.position() = newPosition;
      particle.velocity() = newVelocity;
    }
  }
  
  if (scenario == SCENARIO_FAUCET && PARTICLE::count < MAX_PARTICLES) {
    generateFaucetParticleSet();
  }
  
  updateGrid();
  
  iteration++;
}

///////////////////////////////////////////////////////////////////////////////
// Verlet integration
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::stepVerletBrute(double dt)
{
  if (PARTICLE::count < MAX_PARTICLES && iteration % 12 == 0) {
    generateFaucetParticleSet();
  }
  
  calculateAccelerationBrute();
  
  for (unsigned int i = 0; i < PARTICLE::count; i++) {
    
    PARTICLE& particle = _particles[i];
    
    VEC3D newPosition = particle.position() + particle.velocity()*dt + particle.acceleration()*dt*dt;
    VEC3D newVelocity = (newPosition - particle.position()) / dt;
    
    particle.position() = newPosition;
    particle.velocity() = newVelocity;
  }
  
  iteration++;
  
}

/*
 Calculate the acceleration of each particle using a grid optimized approach.
 For each particle, only particles in the same grid cell and the (26) neighboring grid cells must be considered,
 since any particle beyond a grid cell distance away contributes no force.
*/
void PARTICLE_SYSTEM::calculateAcceleration() {
  
  ///////////////////
  // STEP 1: UPDATE DENSITY & PRESSURE OF EACH PARTICLE
  
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
                  
                  VEC3D diffPosition = particle.position() - neighborParticle.position();
                  
                  double radiusSquared = diffPosition.dot(diffPosition);
                  
                  if (radiusSquared <= h*h)
                    particle.density() += Wpoly6(radiusSquared);
                  
                }
              }
            }
          }
          
          particle.density() *= PARTICLE_MASS;
                        
          // p = k(density - density_rest)
          
          particle.pressure() = GAS_STIFFNESS * (particle.density() - REST_DENSITY);           
        }
      }
    }
  }
  
  ///////////////////
  // STEP 2: COMPUTE FORCES FOR ALL PARTICLES
  
  for (int x = 0; x < (*grid).xRes(); x++) {
    for (int y = 0; y < (*grid).yRes(); y++) {
      for (int z = 0; z < (*grid).zRes(); z++) {
        
        vector<PARTICLE>& particles = (*grid)(x,y,z);
        
        for (int p = 0; p < particles.size(); p++) {
          
          PARTICLE& particle = particles[p];
          
          //cout << "particle id: " << particle.id() << endl;
          
          VEC3D f_pressure, 
          f_viscosity, 
          f_surface, 
          f_gravity = gravityVector * particle.density(),
          colorFieldNormal;
          
          double colorFieldLaplacian;
                    
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
                  
                  //if (particle.id() == neighbor.id()) continue; // SKIPPING COMPARISON OF THE SAME PARTICLE
                  
                  VEC3D diffPosition = particle.position() - neighbor.position();
                  double radiusSquared = diffPosition.dot(diffPosition);
                  
                  if (radiusSquared <= h*h) {
                    
                    VEC3D poly6Gradient, spikyGradient;
                    
                    Wpoly6Gradient(diffPosition, radiusSquared, poly6Gradient);
                    
                    WspikyGradient(diffPosition, radiusSquared, spikyGradient);
                    
                    if (particle.id() != neighbor.id()) {
                      
                      f_pressure += (particle.pressure()/pow(particle.density(),2)+neighbor.pressure()/pow(neighbor.density(),2))*spikyGradient;
                      
                      f_viscosity += (neighbor.velocity() - particle.velocity()) * WviscosityLaplacian(radiusSquared) / neighbor.density();

                    }
                    
                                        
                    colorFieldNormal += poly6Gradient / neighbor.density();
                    
                    colorFieldLaplacian += Wpoly6Laplacian(radiusSquared) / neighbor.density();
                    
                  }
                }
              }
            }
          } // end of neighbor grid cell iteration
          
          f_pressure *= -PARTICLE_MASS * particle.density();
          
          f_viscosity *= VISCOSITY * PARTICLE_MASS;
          
          colorFieldNormal *= PARTICLE_MASS;
          
          
          particle.normal = -1.0 * colorFieldNormal;
          
          colorFieldLaplacian *= PARTICLE_MASS;
          
          
          // surface tension force
          
          double colorFieldNormalMagnitude = colorFieldNormal.magnitude();
          
          if (colorFieldNormalMagnitude > SURFACE_THRESHOLD) {
            
            particle.flag() = true;
            f_surface = -SURFACE_TENSION * colorFieldNormal / colorFieldNormalMagnitude * colorFieldLaplacian;
            
          }
          else {
            particle.flag() = false;
          }
          
          // ADD IN SPH FORCES
          
          particle.acceleration() = (f_pressure + f_viscosity + f_surface + f_gravity) / particle.density();
          
          // EXTERNAL FORCES HERE (USER INTERACTION, SWIRL)
          
          VEC3D f_collision; 
          collisionForce(particle, f_collision);
          
        } 
      }
    }
  }
  
}


void PARTICLE_SYSTEM::collisionForce(PARTICLE& particle, VEC3D& f_collision) {
    
  for (unsigned int i = 0; i < _walls.size(); i++) {
    
    WALL& wall = _walls[i];
    
    double d = (wall.point() - particle.position()).dot(wall.normal()) + 0.01; // particle radius
    
    if (d > 0.0) {
      // This is an alernate way of calculating collisions of particles against walls, but produces some jitter at boundaries
      //particle.position() += d * wall.normal();
      //particle.velocity() -= particle.velocity().dot(wall.normal()) * 1.9 * wall.normal();
      
      particle.acceleration() += WALL_K * wall.normal() * d;
      particle.acceleration() += WALL_DAMPING * particle.velocity().dot(wall.normal()) * wall.normal();
    }
  }
}


double PARTICLE_SYSTEM::Wpoly6(double radiusSquared) { 
    
  static double coefficient = 315.0/(64.0*M_PI*pow(h,9));
  static double hSquared = h*h;
  
  return coefficient * pow(hSquared-radiusSquared, 3);
}

void PARTICLE_SYSTEM::Wpoly6Gradient(VEC3D& diffPosition, double radiusSquared, VEC3D& gradient) {  
    
  static double coefficient = -945.0/(32.0*M_PI*pow(h,9));
  static double hSquared = h*h;
    
  gradient = coefficient * pow(hSquared-radiusSquared, 2) * diffPosition;
}

double PARTICLE_SYSTEM::Wpoly6Laplacian(double radiusSquared) {
    
  static double coefficient = -945.0/(32.0*M_PI*pow(h,9));
  static double hSquared = h*h;
    
  return coefficient * (hSquared-radiusSquared) * (3.0*hSquared - 7.0*radiusSquared);
}

void PARTICLE_SYSTEM::WspikyGradient(VEC3D& diffPosition, double radiusSquared, VEC3D& gradient) {  // 
  
  static double coefficient = -45.0/(M_PI*pow(h,6));
  
  double radius = sqrt(radiusSquared);
  
  gradient = coefficient * pow(h-radius, 2) * diffPosition / radius;
}
 

double PARTICLE_SYSTEM::WviscosityLaplacian(double radiusSquared) {  
  
  static double coefficient = 45.0/(M_PI*pow(h,6));
  
  double radius = sqrt(radiusSquared);
  
  return coefficient * (h - radius);    
}

/*
 Calculate the acceleration of every particle in a brute force manner (n^2).
 Used for debugging.
*/
void PARTICLE_SYSTEM::calculateAccelerationBrute() {
        
  ///////////////////
  // STEP 1: UPDATE DENSITY & PRESSURE OF EACH PARTICLE
  
  for (int i = 0; i < PARTICLE::count; i++) {
    
    // grab ith particle reference
    
    PARTICLE& particle = _particles[i];
    
    // now iteratate through neighbors
    
    particle.density() = 0.0;
            
    for (int j = 0; j < PARTICLE::count; j++) {
      
      PARTICLE& neighborParticle = _particles[j];
      
      VEC3D diffPosition = particle.position() - neighborParticle.position();
      
      double radiusSquared = diffPosition.dot(diffPosition);
      
      if (radiusSquared <= h*h)
        particle.density() += Wpoly6(radiusSquared);
    }
    
    particle.density() *= PARTICLE_MASS;
    
    // p = k(density - density_rest)
    
    particle.pressure() = GAS_STIFFNESS * (particle.density() - REST_DENSITY);
  
    //totalDensity += particle.density();
  }
    
  ///////////////////
  // STEP 2: COMPUTE FORCES FOR ALL PARTICLES
  
  for (int i = 0; i < PARTICLE::count; i++) {
    
    PARTICLE& particle = _particles[i];
    
    //cout << "particle id: " << particle.id() << endl;
    
    VEC3D f_pressure, 
    f_viscosity, 
    f_surface, 
    f_gravity(0.0, particle.density() * -9.80665, 0.0),
    n, 
    colorFieldNormal,
    colorFieldLaplacian; // n is gradient of colorfield
    //double n_mag;
    
    for (int j = 0; j < PARTICLE::count; j++) {
      PARTICLE& neighbor = _particles[j];
      
      VEC3D diffPosition = particle.position() - neighbor.position();
      VEC3D diffPositionNormalized = diffPosition.normal(); // need?
      double radiusSquared = diffPosition.dot(diffPosition);
      
      if (radiusSquared <= h*h) {
        
                
        if (radiusSquared > 0.0) {
          
          //neighborsVisited++;
          //cout << neighborsVisited << endl;
          
          //cout << neighbor.id() << endl;
          
          VEC3D gradient;
                    
          Wpoly6Gradient(diffPosition, radiusSquared, gradient);
                    
          f_pressure += (particle.pressure() + neighbor.pressure()) / (2.0 * neighbor.density()) * gradient;
          
          colorFieldNormal += gradient / neighbor.density();
        }
        
        f_viscosity += (neighbor.velocity() - particle.velocity()) * WviscosityLaplacian(radiusSquared) / neighbor.density();
        
        colorFieldLaplacian += Wpoly6Laplacian(radiusSquared) / neighbor.density();
      }
      
    }
    
    f_pressure *= -PARTICLE_MASS;
    
    //totalPressure += f_pressure;
    
    f_viscosity *= VISCOSITY * PARTICLE_MASS;
    
    colorFieldNormal *= PARTICLE_MASS;
    
    particle.normal = -1.0 * colorFieldNormal;
    
    colorFieldLaplacian *= PARTICLE_MASS;
    
    // surface tension force
    
    double colorFieldNormalMagnitude = colorFieldNormal.magnitude();
    
    if (colorFieldNormalMagnitude > surfaceThreshold) {
      
      particle.flag() = true;
      f_surface = -SURFACE_TENSION * colorFieldLaplacian * colorFieldNormal / colorFieldNormalMagnitude;
      
    }
    else {
      particle.flag() = false;
    }
    
    // ADD IN SPH FORCES
    
    particle.acceleration() = (f_pressure + f_viscosity + f_surface + f_gravity) / particle.density();
    
    
    // EXTERNAL FORCES HERE (USER INTERACTION, SWIRL)
    
    VEC3D f_collision;
    
    collisionForce(particle, f_collision);    

  }
}

