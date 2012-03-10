#include "PARTICLE_SYSTEM.h"


///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
PARTICLE_SYSTEM::PARTICLE_SYSTEM() :
  grid(10,10,10)
{
  f_gravity = VEC3F(0, -9.8, 0);
  
  
  //std::cout << grid.size() << "\n";
  //gridKey key = std::tr1::make_tuple(2,3,4);
  //vector<PARTICLE> *hi = &grid[key];
  //if (gridKey.find("
   
  /* // 3d fluid box in middle
  for (int x = -3; x < 4; x++) {
    for (int y = 9; y < 12; y++) {
      for (int z = -3; z < 4; z++) {
        _particles.push_back(PARTICLE(VEC3F(x*0.02,y*0.02,z*0.02)));
      }
      
    }
  }
   
   */
  
  /*
  _particles.push_back(PARTICLE(VEC3F(0, 0, 0)));
  _particles.push_back(PARTICLE(VEC3F(0, 0.02, 0)));
  _particles.push_back(PARTICLE(VEC3F(0, 0.04, 0)));
*/

  bool offset = false;
  
  for (float y = 0; y < BOX_SIZE/2.0; y+= 0.02) {
    for (float x = BOX_SIZE/4.0; x < BOX_SIZE/2.0; x += 0.02) {
      for (float z = -BOX_SIZE/2.0; z < BOX_SIZE/2.0; z+= 0.02) {
        _particles.push_back(PARTICLE(VEC3F((offset ? -0.01 : 0) + x, y,z)));
      }
    }
    offset = !offset;
  }
   
  
  printf("simulating %d particles\n", (int)_particles.size());
  
  _walls.push_back(WALL(VEC3F(0,0,1), VEC3F(0,0,-BOX_SIZE/2.0)));
  _walls.push_back(WALL(VEC3F(1,0,0), VEC3F(-BOX_SIZE/2.0,0,0)));
  _walls.push_back(WALL(VEC3F(-1,0,0), VEC3F(BOX_SIZE/2.0,0,0)));
  _walls.push_back(WALL(VEC3F(0,1,0), VEC3F(0,-BOX_SIZE/2.0,0)));
  _walls.push_back(WALL(VEC3F(0,0,-1), VEC3F(0,0,BOX_SIZE/2.0)));
  //_walls.push_back(WALL(VEC3F(0,-1,0), VEC3F(0,BOX_SIZE/2.0,0)));

  
  /*
  
  std::cout << grid(3,4,5).size() << "\n";
  
  PARTICLE hi = PARTICLE(VEC3F(33,44,52));
  grid(3,4,5).push_back(hi);
  _particles.push_back(hi);
  
  PARTICLE hi1 = PARTICLE(VEC3F(34,44,52));
  grid(3,4,5).push_back(hi1);
  _particles.push_back(hi1);
  
  std::cout << grid(3,4,5).size() << "\n";
  
  
  std::cout << grid.size() << "\n";
  //calculateAcceleration();
   
   
   
   */
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
  for (unsigned int x = 0; x < _particles.size(); x++)
    _particles[x].draw();

  
  glutWireCube(BOX_SIZE);
  
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
  unsigned int numParticles = _particles.size();
  
  //VEC3F systemAcceleration[numParticles];
  
  calculateAcceleration();
  
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


void PARTICLE_SYSTEM::calculateAcceleration() {
  
  unsigned int numParticles = _particles.size();
  
  // to conaint a vector of neighbors belonging to each particle (same index as _particles)
  vector<PARTICLE> **neighborArray = new vector<PARTICLE>*[numParticles]; // make this vector instead?
  
  
  ///////////////////
  // STEP 1: UPDATE DENSITY & PRESSURE OF EACH PARTICLE
  
  for (int i = 0; i < numParticles; i++) {
    
    // grab ith particle reference
    
    PARTICLE& particle = _particles[i];
    
    //vector<PARTICLE> totalNeighborParticles;
    
    
    /*
    
    // get grid position from position (hash it basically)
    
    VEC3F gridPosition = particle.position() / h;
    
    // get the neighbors of this particle
    
    getNeighborParticles(totalNeighborParticles, (int)gridPosition.x, (int)gridPosition.y, (int)gridPosition.z);
    
    
    */
    
    //totalNeighborParticles = _particles;
     
    //cout << neighborArray[i]->size() << "\n";
    
    neighborArray[i] = &_particles;
    
    // now iteratate through neighbors
    
    unsigned int numNeighborParticles = neighborArray[i]->size();
    
    float density = 0.0;
    
    for (int j = 0; j < numNeighborParticles; j++) {
      
      PARTICLE& neighborParticle = (*neighborArray[i])[j];
      
      // 1st get density
      
      // r squared is x^2 + y^2 + z^2
      
      VEC3F diffPosition = particle.position() - neighborParticle.position();
      
      density += PARTICLE_MASS * Wpoly6(diffPosition);
    }
    
    particle.density() = density;
    
    // then from density get the pressure    p = k(density - density_rest)
    
    particle.pressure() = PARTICLE_K * (particle.density() - PARTICLE_REST_DENSITY);
    
    //cout << particle.pressure() << endl;
    
  }
  
  
  ///////////////////
  // STEP 2: COMPUTE FORCE VECTORS FOR ALL PARTICLES
  
  for (int i = 0; i < numParticles; i++) {
    
    VEC3F f_pressure, f_viscosity, f_surface, n; // n is gradient of colorfield
    float n_mag;
    PARTICLE& particle = _particles[i];
    vector<PARTICLE>& neighbors = *(neighborArray[i]);
        
    unsigned int numNeighborParticles = neighbors.size();
    
    //cout << "numNeighborParticles:" << numNeighborParticles << endl;
    
    //cout << _particles[0].density() << " " << _particles[1].density() << endl;
    //cout << neighbors[0].density() << " " << neighbors[1].density() << endl;
    //_particles[0].density() = 4.0;
    //cout << _particles[0].density() << " " << _particles[1].density() << endl;
    //cout << neighbors[0].density() << " " << neighbors[1].density() << endl;

    for (int j = 0; j < numNeighborParticles; j++) {
      
      PARTICLE& neighbor = neighbors[j];
      //cout << "neighbor position: " <<  neighbor.position() << "\n";
      //cout << "particle position: " << particle.position() << "\n";
      
      VEC3F diffPosition = particle.position() - neighbor.position();
      VEC3F diffPositionNormalized = diffPosition.normal();
      
      
      //cout << "diff position: " << diffPosition << endl;
      //cout << "particle pressure:" << particle.pressure() << " density:" << particle.density() << endl;
      //cout << "neighbor pressure:" << neighbor.pressure() << " density:" << neighbor.density() << endl;

      //cout << "particle pos: " << particle.position() << endl;
      //cout << "neighbor pos: " << neighbor.position() << endl;
      //cout << "diffPosition: " << diffPosition << endl;
      
      f_pressure -= PARTICLE_MASS * (particle.pressure() + neighbor.pressure()) / (2 * neighbor.density()) * Wspiky_1(diffPosition) * diffPositionNormalized;
      
      //cout << f_pressure << "\n";
      
      f_viscosity += PARTICLE_MASS * (neighbor.velocity() - particle.velocity()) / neighbor.density() * Wviscosity_2(diffPosition) * diffPositionNormalized;
      
      n += PARTICLE_MASS / neighbor.density() * Wpoly6_1(diffPosition) * diffPositionNormalized;
      
    }
    
    f_viscosity *= PARTICLE_MEW;
    
    // check if |n| > threshold...if so continue finding surface force
    
    n_mag = n.magnitude();
    
    //cout << "n_mag: " << n_mag << endl;
    
    //cout << n_mag << endl;
    
    if (n_mag > N_SURFACE_THRESHOLD) {  /// FALSE FLASE FALSE FALSE 
      
      VEC3F colorField_2; // laplacian of colorfield
      
      for (int j = 0; j < numNeighborParticles; j++) {
        
        PARTICLE& neighbor = neighbors[j];
        
        VEC3F diffPosition = particle.position() - neighbor.position();
      
        colorField_2 += PARTICLE_MASS / neighbor.density() * Wpoly6_2(diffPosition) * diffPosition.normal();
        
      }
      
      particle.flag() = true;
      
      f_surface = -PARTICLE_THETA * colorField_2 * (n / n_mag);
      
      //cout << "f_surface: " << f_surface << endl;
      
    }
    else {
      particle.flag() = false;
    }
    
    // ADD IN SPH FORCES
    
    particle.acceleration() = (f_pressure + f_viscosity + f_surface) / particle.density();
    
    
    // EXTERNAL FORCES HERE (USER INTERACTION, SWIRL)
    
    VEC3F f_collision;
    
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
    
    
    
    particle.acceleration() += (f_gravity + f_collision) / PARTICLE_MASS;
  }
  
  delete[] neighborArray;
}

void PARTICLE_SYSTEM::getNeighborParticles(vector<PARTICLE>& totalNeighborParticles, int x, int y, int z) {
  
  vector<PARTICLE>& neighborParticles = grid(x, y, z);
  totalNeighborParticles.insert(totalNeighborParticles.end(), neighborParticles.begin(), neighborParticles.end());
  
  if (x > 0) {
    vector<PARTICLE>& neighborParticles = grid(x - 1, y, z);
    totalNeighborParticles.insert(totalNeighborParticles.end(), neighborParticles.begin(), neighborParticles.end()); 
  }
  
  if (x < grid.xRes() - 1) {
    vector<PARTICLE>& neighborParticles = grid(x + 1, y, z);
    totalNeighborParticles.insert(totalNeighborParticles.end(), neighborParticles.begin(), neighborParticles.end()); 
  }
  
  if (y > 0) {
    vector<PARTICLE>& neighborParticles = grid(x, y - 1, z);
    totalNeighborParticles.insert(totalNeighborParticles.end(), neighborParticles.begin(), neighborParticles.end()); 
  }
  
  if (y < grid.yRes() - 1) {
    vector<PARTICLE>& neighborParticles = grid(x, y + 1, z);
    totalNeighborParticles.insert(totalNeighborParticles.end(), neighborParticles.begin(), neighborParticles.end()); 
  }
  
  if (z > 0) {
    vector<PARTICLE>& neighborParticles = grid(x, y, z - 1);
    totalNeighborParticles.insert(totalNeighborParticles.end(), neighborParticles.begin(), neighborParticles.end()); 
  }
  
  if (z < grid.zRes() - 1) {
    vector<PARTICLE>& neighborParticles = grid(x, y, z + 1);
    totalNeighborParticles.insert(totalNeighborParticles.end(), neighborParticles.begin(), neighborParticles.end()); 
  }
  
}

float PARTICLE_SYSTEM::Wpoly6(VEC3F& r) {  // checked
  
  float radiusSquared = r.dot(r);
  
  if (radiusSquared < 0 || radiusSquared > pow(h,2)) {  // IMPORTANT! since we are checking r^2, need to check against h^2
    return 0.0;
  }
  else {
    return 315.0/(64*pi*pow(h,9)) * pow( pow(h,2)-radiusSquared, 3);
  }
  
}

float PARTICLE_SYSTEM::Wpoly6_1(VEC3F& r) {  // checked
  
  float radius = r.magnitude();
  
  if (radius < 0 || radius > h) {
    return 0.0;
  }
  else {
    return 315.0/(64*pi*pow(h,9)) * -6.0*radius*pow(pow(h,2) - pow(radius,2), 2); // optimize here
  }
}

float PARTICLE_SYSTEM::Wpoly6_2(VEC3F& r) { // checked
  
  float radiusSquared = r.dot(r);
  
  if (radiusSquared < 0 || radiusSquared > pow(h,2)) {  // IMPORTANT! since we are checking r^2, need to check against h^2
    return 0.0;
  }
  else {
    return 315.0/(64*pi*pow(h,9)) * -6.0*(pow(h,2)-5.0*radiusSquared)*(pow(h,2)-radiusSquared) - 12.0*pow(pow(h,2)-radiusSquared, 2);
  }
  
}

float PARTICLE_SYSTEM::Wspiky_1(VEC3F& r) {  // checked
  
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

float PARTICLE_SYSTEM::Wviscosity_2(VEC3F& r) {  // checked
  
  float radius = r.magnitude();
  
  if (radius < 0 || radius > h) {
    return 0.0;
  }
  else {
    return 45/(pi*pow(h,6))*(h - radius);  // fixed to match the paper's laplacian of viscosity
  }
  
}
