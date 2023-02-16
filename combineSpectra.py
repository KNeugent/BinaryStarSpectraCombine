import numpy as np
from scipy.interpolate import interp1d


def convert_vac_to_air(vac_wave):
    """
    Converts wavelength from vacuum wavelength to air
    Conversion comes from Morton (1991, ApJS, 77, 119)

        Parameters:
            vac_wave ([float]): array of wavelength values in vacuum

        Returns:
            air_wave ([float]): array of wavelength values in air
    """
    air_wave = vac_wave / (
        1.0 + 2.735182e-4 + 131.4182 / (vac_wave**2) + 2.76249e8 / (vac_wave**4)
    )

    return air_wave


def calc_radius(lum, temp):
    """
    Calculates the radius given the temperature and luminosity of a star

        Parameters:
            lum (float): luminosity of star (ergs)
            temp (float): temperature of star (K)

        Returns:
            radius (float): radius of star (cm)
    """
    # Stefan-Boltzmann constant
    sig = 5.6e-5  # erg cm^-2 K^-4 K^-1
    radius = np.sqrt(lum / (4 * np.pi * temp**4**sig))

    return radius


def calc_lum(Mv, BC):
    """
    Calculates the luminosity given an M_V and bolometric correction.

        Parameters:
            Mv (float): M_V / absolute magnitude
            BC (float): Bolometric Correction

        Returns:
            lum (float): Luminosity (ergs)
    """
    lum = 10 ** (((Mv + BC) - 4.7) / -2.5) * (3.85e33)

    return lum


def interp_spectra(wavelength, flux, n_points):
    """
    Smoothes / interpolates the spectrum to be the requested number of points.

        Parameters:
            wavelength ([float]): wavelength array
            flux ([float]): float array
            n_points (int): total number of points requested in interpolated array

        Returns:
            wavelength_interp ([float]): interpolated wavelength array
            flux_interp ([float]): interpolated flux array
    """
    # create a linear interpolation for wavelength and flux
    lin_interp = interp1d(wavelength, flux)
    # create a linear spacing between the max and min wavelenth with
    # the requested number of points
    wavelength_interp = np.linspace(min(wavelength), max(flux), n_points)
    # interpolate the flux
    flux_interp = lin_interp(wavelength_interp)

    return wavelength_interp, flux_interp


def redden_flux(A_v, wavelength, flux):
    """
    Reddens the flux given the input A_v value.

        Parameters:
            A_v (float): extinction measured in the V band
            wavelength ([float]): wavelength array
            flux ([float]): flux array

        Returns:
            reddened_flux ([float]): reddened flux array
    """
    # assumes normal Cardelli et al. (1989) reddening
    R_v = 3.1

    # transformations below come from Carpenter (2001)
    x = 1 / (wavelength * 10 ** (-4))
    y = x - 1.82

    a_x = (
        1
        + 0.17699 * y
        - 0.50447 * y**2
        - 0.02427 * y**3
        + 0.72085 * y**4
        + 0.01979 * y**5
        - 0.77530 * y**6
        + 0.32999 * y**7
    )
    b_x = (
        1.41338 * y
        + 2.28305 * y**2
        + 1.07233 * y**3
        - 5.38434 * y**4
        - 0.62251 * y**5
        + 5.30260 * y**6
        - 2.09002 * y**7
    )

    A_lambda = (a_x + (b_x / R_v)) * A_v
    z = 10 ** (A_lambda / 2.5)
    reddened_flux = flux / z

    return reddened_flux


def scale_luminosity(luminosity, flux, multiplier):
    """
    Scales the luminosity by a multiplier to account for an increase
    or decrease in mass.

        Parameters:
            luminosity (float): luminosity (ergs)
            flux ([float]): flux array
            multiplier (float): luminosity multiplier (logL)

        Returns:
            scaled_flux ([float]): flux array given the new luminosity
    """
    # convert luminosity to logL
    logL = np.log10(luminosity)
    # determine the difference
    scale_1 = logL - multiplier
    # take out of log units
    scale_2 = 10.0**scale_1
    scaled_flux = flux / scale_2

    return scaled_flux


def combine_spectra(star1_flux, star2_flux):
    """
    Combines the flux of two spectra by adding them together.

        Parameters:
            star1_flux ([float]): flux array of star 1
            star2_flux ([float]): flux array of star 2

        Returns:
            combined_flux ([float]): combined flux of star 1 and star2
    """
    combined_flux = star1_flux + star2_flux

    return combined_flux


# The remaining helper functions are specific to the MARCS and BSTAR06 models


def rad_lum_from_header(filename):
    """
    Extract the radius and luminosity from the MARCS spectra header

        Parameters:
            filename (string): file name of MARCS file

        Returns:
            radius (float): radius of MARCS model star
            luminosity (float): luminosity of MARCS model star
    """
    ## RSG header to determine radius
    with open(filename + ".mod") as file:
        for i, line in enumerate(file):
            if i == 7:
                radius = float(line[0:12])
            if i == 8:
                luminosity = float(line[0:12])
                break

    return radius, luminosity


def extract_BSTAR_vals(filename):
    """
    Extract the relevant information from the BSTARO6 model header
    and use this information to convert the flux to the same
    units as the MARCS models.

        Parameters:
            filename (string): name of BSTAR06 model file

        Returns:
            wavelength ([float]): wavelength array
            flux ([float]): flux array
    """
    # reads the first line of the file
    with open(filename, "r") as f:
        inline = f.readline().strip("\n")

    # saves the relevant information into variables
    cl = line[1]
    temperature = float(line[3:8])
    Mv = float(line[25:29])
    BC = float(line[33:37])
    g = line[41:44]

    # determines the class of star based on the number
    # in the first line
    class_dict = {"5": "V", "1": "I", "3": "III"}
    ty = class_dict[cl]

    # calculate the luminosity given Mv and BC
    luminosity = calc_lum(Mv, BC)
    # calcualte the radius given the luminosity
    radius = calc_radius(luminosity, temperature)

    # read in the actual model
    BSTAR_model = np.loadtxt(filename, comments="#")
    # interpolate the spectrum so it can later be joined with the MARCS model
    # also, multiply the flux by the radius to convert the flux to the
    # same units as the MARCS models
    wavelength, flux = interp_spectra(BSTAR_model[:, 0], BSTAR_model[:, 1] * radius, 7000)

    return wavelength, flux


def extract_MARCS_vals(filename):
    """
    Extract the relevant information from the MARCS model header
    and get spectrum ready to be combined with BSTAR06 models.

        Parameters:
            filename (string): name of MARCS model file

        Returns:
            wavelength ([float]): wavelength array
            flux ([float]): flux array
            radius (float): radius from header
            luminosity (float): luminosity from header
    """
    file_flux = filename + ".flx"
    file_mod = filename + ".mod"

    # read in the wavelength values for the MARCS models
    MARCS_wave = np.loadtxt("MARCS/flx_wavelengths.vac")
    MARCS_flux = np.loadtxt(file_flux)
    # convert from vacuum to air
    MARCS_waveAir = convert_vac_to_air(MARCS_wave)

    # shorten MARCSarray to only include visible wavelength range
    wave_sm = MARCS_waveAir[18017:40766]
    flux_sm = MARCS_flux[18017:40766] * R_R

    # interpolate
    wavelength, flux = interp_spectra(wave_sm, flux_sm, 7000)

    # get radius and luminosity from header
    # to later be used to scale the flux appropriately
    radius, luminosity = rad_lum_from_header(file_mod)

    return wavelength, flux, radius, luminosity


def redden_scale(filename, lum_multiplier, reddening):
    """
    Helper function to redden and scale the luminosity / flux values
    after reading in the spectrum.

        Parameters:
            filename (string): name of MARCS file
            lum_multiplier (float): luminosity multiplier (logL)
            reddening (float): A_v value to redden the spectrum by

        Returns:
            wavelength ([float]): wavelength array
            flux_red ([float]): reddened and scaled flux array
    """
    wavelength, flux, radius, luminosity = extract_MARCSvals(filename)
    flux_scaled = scale_luminosity(luminosity, flux, lum_multiplier)
    flux_red = redden_flux(reddening, wavelength, flux_scaled)

    return wavelength, flux_red


def main():
    # the MARCS models are created for unreddened, 15Mo stars
    # this code will generate red supergiants of different masses
    # by scaling the flux based on the luminosity and also apply
    # a reddening correction

    # the following luminosity multiplier corresponds to the
    # following masses: 10Mo, 15Mo, 20Mo, 25Mo
    lum_multiplier = [4.2, 1, 5.2, 5.4]

    # A_v / reddening value
    # appropriate values range from 0 - 2.8 in 0.2 increments
    reddening = 1.2

    # STAR 1
    # in this example, this is a red supergiant from MARCS models
    star1_filename = "s4000_g+0.0_m15._t05"
    star1_wave, star1_flux = redden_scale(star1_filename, lum_multiplier, reddening)

    # STAR 2
    # in this example, this is a B star from the BSTAR06 models
    star2_filename = "BG20000g400_150.11"
    star2_wave, star2_flux = extract_BSTAR_vals(star2_filename)

    flux_combined = combine_spectra(star1_flux, star2_flux)

    outputFileName = "RSG_Bstar_combinedSpec.txt"

    # save the output spectrum
    np.savetxt(outputFileName, np.array(list(zip(star2_wave, flux_combined))))


if __name__ == "__main__":
    main()
