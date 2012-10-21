For our final project in Theodore Kim's Physics-Based Simulation class we were asked to implement a large, non-trivial technique. My project is based on Matthias Müller's (et al.) research of applying Smoothed Particle Hydrodynamics (SPH) to fluid simulation - essentially simulating fluid as thousands of small particles, rather than the traditional Eularian grid-based approach.

Feel free to email me and I'll tell you more about this project!

Basic Controls:
Click & drag to rotate the scene
a - Toggle the physics animation on/off (upon starting animation is OFF)
1 - Load the "Dam Break" scenario (initial scenario)
2 - Load the "Cube" scenario
3 - Load the "Faucet" scenario
t - Toggle tumble mode on/off (the water will react to the tilt of the camera)

Advanced Controls:
g - Toggle grid visible (grid used for neighbor particle optimization)
/ - Toggle gravity on/off
s - Toggle coloring of surface purple (particles on the surface have a surface tension force applied)
f - Output current FPS to stdout

Here is more information on the paper I used:

Particle-Based Fluid Simulation for Interactive Applications
Matthias Müller, David Charypar and Markus Gross
http://www.matthiasmueller.info/publications/sca03.pdf