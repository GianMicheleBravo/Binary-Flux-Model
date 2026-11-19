# Binary-Flux-Model
A model to compute the flux of an eclipsing binary system at any time, at any wavelength. Optimized for speed and for a minimum number of free parameters, and suitable for detached systems without significant ellipsoidal stretching. Depending on the setup, this model can calculate the flux of 10⁴ - 10⁵ data points per second.

## Quick start guide:
- Clone the repository to your machine, in any directory accessible by your scripts.
- Add the code ```import binary_flux_model``` into your script.
- Optional step: use the code ```binary_flux_model.plot_eclipsing_binary_flux(2021, 4, 2024, 80, 1e-8, 0.2, 0.3, 15000, 12000)``` to produce the following plot:
<img width="500" height="475" alt="demo_plot_2" src="https://github.com/user-attachments/assets/164e0cd4-8c6b-49dc-9f27-dc1be505fc14" />

- Optional step: watch the [demo animation on PCloud](https://e.pcloud.link/publink/show?code=XZ0UDV7ZUgMOIOyACMmB2UfngsbkUYGLLT3V#) which uses frames produced by the ```plot_eclipsing_binary_flux``` function.
- Get the flux of using the ```binary_flux_model.eclipsing_binary_flux``` 

## Getting the flux:
The flux of an eclipsing binary system is obtained as following:
```
binary_flux_model.eclipsing_binary_flux(t, lambd, P, t0, i_deg, theta_max, f1, f2, T1, T2, spectrum_function_1=Planck_law, spectrum_function_2=Planck_law):
    """Return the light emitted by an eclipsing binary. 
        This model assumes the objects as uniform, opaque spheres in a circular orbit.
    Input:
        t: time, array-like or scalar
        lambd: wavelength, array-like with or scalar, larger than 0, default in nm
        P: period, array-like with or scalar, larger than 0, same unit as t
        t0: time of primary eclipse, array-like with or scalar, same unit as t
        theta_max: maximum projected separation angle of the binary, array-like with or scalar,
            default in radians
        i_deg: angle between the plane vector of the system and the line of sight,
            array-like or scalar, between 0 and 90, must be in degrees
        f1: the ratio between the primary radius and the binary distance, array-like or scalar,
            between 0 and 1
        f2: the ratio between the secondary radius and the binary distance,array-like or scalar,
            between 0 and 1
        T1: temperature of the primary, array-like or scalar, default in K
        T2: temperature of the secondary, array-like or scalar, default in K

        spectrum_function_1: the function used to get flux/solid angle of the primary:
            spectrum_function_1(lambd, T) -> flux
            Default: Planck_law:
                lambd_nm, T_K -> flux in Jy/sr

        spectrum_function_2: the function used to get flux/solid angle of the secondary:
            spectrum_function_2(lambd, T) -> flux
            Default: Planck_law:
                lambd_nm, T_K -> flux in Jy/sr

    Output:
        out: the flux of the system, unit specificed by
            spectrum_function_1 and spectrum_function_2, default in Jy

    Example:
        >>> t = [0, 10, 20, 30, 40]
        >>> lambd_nm = 550 # green
        >>> P = 40
        >>> t0 = 10
        >>> i_deg = 85
        >>> theta_max_rad = 1/21600000 # 1 milliarcsecond
        >>> f1 = 0.3
        >>> f2 = 0.4
        >>> T1_K = 6000
        >>> T2_K = 5000
        >>> eclipsing_binary_flux(time, lambd_nm, P, t0, i_deg, theta_max_rad, f1, f2, T1_K, T2_K)
        array([40906.9806891 , 17367.02548793, 40906.9806891 , 31138.02885213,
            40906.9806891 ])"""
```
This is the single most important function of the package, as it gives you the system's flux, that you can compare directly to measurements. For ease of use, any of its input parameters (except for the spectral functions) can be both scalars, lists or arrays, as long as the input have compatible shapes. In simple words, the model works as following:
1. The flux/solid angle is calculated for both objects in the binary for the input wavelength(s) using ```spectrum_function_1``` and ```spectrum_function_2``` respectively. 
2. The solid angle of both objects as visible from Earth is calculated for the input time(s).
3. The visible solid angle of both objects is multiplied by their respective flux/solid angle to get the total flux of the system.
 
This setup makes it possible for the model to compute both the SED (Spectral Energy Distribution) and light curve of an eclipsing binary system, for any time and wavelength.

### Changing the spectral models and units:
By default, the model uses the following units:
- ```t```: any time unit
- ```lambd```: nanometers
- ```P```: same unit as t
- ```t0```: same unit as t
- ```i_deg```: degrees. This unit cannot be changed
- ```theta_max```: radians
- ```f1```: unitless. This cannot be changed
- ```f2```: unitless. This cannot be changed
- ```T1```: Kelvin
- ```T2```: Kelvin
Using these units and ```binary_flux_model.Planck_law``` as a spectrum function, the flux/solid angle is calculated in Janskys/steradians (Jy/Sr), the solid angle is calculated in steradians, and the output flux is in Janskys. However, if the objects in the binary do not emit according to Planck's law, or we want different input and output units, we can change this function. As an example, let us try this:
```
def new_spectrum_function(lambd_micrometers, T_Celsius):
    """A made up function that return the flux in counts/s/milliarcsecond^2
    from the wavelenght in micrometers and the temperature in Celsius"""

    T_Kelvin = T_Celsius + 273.15
    flux_counts_per_s_per_mas2 = 2**-(lambd_microns - (T_Kelvin / 100))
    flux_counts_per_s_per_mas2 *= 1000 * T_Kelvin
    return flux_counts_per_s_per_mas2
```
Yes, that is a Gaussian function. Yes, it does not make sense. Yes, it uses 100 and 1000 as scaling parameters. Yes, you can make a better spectrum function yourself. The point is that any function works, as long as it takes wavelength as an input, and gives output in flux/solid angle. This works with any unit for wavelength (if you want to be cheeky, you can even sneak frequency measurements through the ```eclipsing_binary_flux``` function and to a custom built spectrum function, and the function will not notice!), temperature, separation angle and flux. We now run the flux model by putting our ```new_spectrum_function``` as its ```spectrum_function_1``` and ```spectrum_function_1``` inputs:
```
binary_flux_model.eclipsing_binary_flux(t, lambd, P, t0, i_deg, theta_max, f1, f2, T1, T2, spectrum_function_1=new_spectral_function, spectrum_function_2=new_spectral_function)
```
We now have an eclipsing star model that uses a new model spectra and has new input and output units. In this case, the input ```lambd``` needs to be in micrometers, the input ```theta_max``` needs to be in milliarcseconds (so that the solid angle calculations become in milliarseconds^2), and the input temperatures needs to be in Kelvin. The output flux is then in counts per second. If desired, ```spectrum_function_1``` and ```spectrum_function_2``` can be different functions, for example if the first object in the binary is a star and the second is a planet.


## Model assumptions:
In order to reduce the number of free parameter, the model makes the following assumptions:
- Both objects are spherical, opaque, and have an homogeneous surface.
- The binary orbit is circular
- Unless any other option is specified, both objects are modeled as blackbody emitters that follow Planck's law.
