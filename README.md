# Combining the Spectra for Binary Stars

This program combines model spectra from two different types of stars to create a final combined binary-star spectrum. It also uses basic stellar astrophysics equations to calculate how the wavelength-dependent flux would change in the following circumstances:
* The mass (luminosity) of a star of a certain temperature and surface gravity is changed.
* A reddening correction is applied.

This program formed the basis of [Neugent et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018AJ....156..225N/abstract) and has been applied to the red supergiant [MARCS models](https://marcs.astro.uu.se) and the B-type star [BSTAR06 models](http://tlusty.oca.eu/Tlusty2002/tlusty-frames-BS06.html). However, it can be generalized to combine model or observed spectra from any two types of stars. The helper functions within the code can also be used to redden a single spectrum, or to transform the wavelength-dependent flux to what you would expect for a more or less massive star of the same temperature and surface gravity.

## Using this code

### Dependencies

The imported packages are `numpy` and `scipy`.

This code has been tested using `python 3.7.3`, `numpy 1.18.2`, and `scipy 1.3.3`.

### Running the code

This code is currently highly specialized to using the MARCS and BSTAR06 models. However, with some edits, it could be used to combine model or observed spectra from any set of stars. 

To run the code on MARCS and BSTAR06 models, the only values that need to be changed are located in the main method. Here you can change:
* the luminosity multiplier (used to change the mass of the red supergiant)
* the reddening (A_v)
* the number of points to interpolate the spectra over
* the input filenames for the two stars (currently defined as one MARCS model and one BSTAR06 model)
* the names of the output spectra (saved as txt files with two columns: wavelength and flux)

Example input files are provided as well. The MARCS model for a 4000K, log g = 0.0, 15Mo, v_turb = 5 km/s star is provided:
* [s4000_g+0.0_m15._t05.flx](https://github.com/KNeugent/BinaryStarSpectraCombine/blob/main/s4000_g%2B0.0_m15._t05.flx) - MARCS model file of flux values
* [s4000_g+0.0_m15._t05.mod](https://github.com/KNeugent/BinaryStarSpectraCombine/blob/main/s4000_g%2B0.0_m15._t05.mod) - MARCS model file of model values (physical properties, etc.)
* [flx_wavelengths.vac](https://github.com/KNeugent/BinaryStarSpectraCombine/blob/main/flx_wavelengths.vac) - MARCS model wavelengths in vacuum (a helper function is in the code to convert wavelengths to air)

The BSTAR06 model for a 20000K, log g = 4.0, v_rot = 150 km/s star is provided:
* [BG20000g400_150.11](https://github.com/KNeugent/BinaryStarSpectraCombine/blob/main/BG20000g400_150.11) - BSTAR06 model file with wavelength, flux, and physical properties in commented first line.
* [BStars.txt](https://github.com/KNeugent/BinaryStarSpectraCombine/blob/main/BStars.txt) - A list of typical B-star temperatures, spectral types, absolute magnitudes (M_v), and bolometric corrections that can be used to edit the input flux to apply to either a dwarf, giant, or supergiant star.

### Outputted files

The program outputs three txt files ([RSG_single_spectrum.txt](https://github.com/KNeugent/BinaryStarSpectraCombine/blob/main/RSG_single_spectrum.txt), [Bstar_single_spectrum](https://github.com/KNeugent/BinaryStarSpectraCombine/blob/main/Bstar_single_spectrum.txt), and [RSG_B_combined_spectrum.txt](https://github.com/KNeugent/BinaryStarSpectraCombine/blob/main/RSG_B_combined_spectrum.txt)) that have the wavelength and flux spanning from 3200 Angstroms - 1 micron of the transformed models.

There is additionally a very small helper program, [plotSpectra.py](https://github.com/KNeugent/BinaryStarSpectraCombine/blob/main/plotSpectra.py) that plots the output files and produces an eps plot, as shown below:

![RSGbinarySpectrum](https://github.com/KNeugent/BinaryStarSpectraCombine/blob/main/RSGbinarySpectrum.jpg)