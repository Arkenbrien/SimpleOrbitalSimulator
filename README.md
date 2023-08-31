# SimpleOrbitalSimulator
Simulates a spherical projectile with a starting velocity above the earth

This is a simple 2D simulator banged out in an afternoon just for fun. Uses gravity and atmospheric drag from a single celestrial body.

I give two programs: planar and spherical. Planar launchs a projectile above a flat surface (with atmosphere), sphrical does the same except with a sphere.

Asumptions of the model:
  - Atmospheric density (affecting drag) was calculated using matlab's atmosisa function, it only goes up to 20000 meters. Program assumes that air doesn't exist above 20000. Fix later? Linear extrapolation when?
  - Earth is simulated as a perfect sphere.
  - I am a decent programer.
  - Weird things happen if the dt AND speed is set too high. If the program does something weird lower the dt to below 1.

Makes a bunch of plots, plus there's an animation if you want.

Allows to vary:
- Radius/mass of cannon ball
- Initial x/z velocities
- Initial altitude above the earth
- Drag coefficient
