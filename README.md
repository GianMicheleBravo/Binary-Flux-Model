# Eclipsing-binary-model
A model to compute the flux of an eclipsing binary system at any time, at any wavelength. Optimized for speed and for a minimum number of free parameters, and suitable for detached systems without significant ellipsoidal stretching. Depending on the setup, this model can calculate the flux of 10⁴ - 10⁵ data points per second.

Installation guide:
- Clone the repository to your machine.
- Add the code ```import binary_flux_model``` into your code. 

Model assumptions:
- Both objects are spherical, opaque, and have an homogeneous surface.
- The binary orbit is circular
- Unless any other option is specified, both objects are modeled as blackbody emitters.
