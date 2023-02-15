import numpy as np
from scipy.interpolate import interp1d

def convert_vac_to_air (vac_wave):
    # convert vaccuum to air
    # needed for MARCS models
    air_wave = vac_wave/(1.0+2.735182E-4+131.4182/(vac_wave**2)+2.76249E8/(vac_wave**4))
    return air_wave

def calc_radius (lum, temp):
    sig = 5.6E-5 # erg cm^-2 K^-4 K^-1
    radius = np.sqrt(lum/(4*np.pi*temp**4**sig))
    return radius

def calc_lum (Mv, BC):
    lum = 10**(((Mv+BC)-4.7)/-2.5)*(3.85E33)
    return lum

def interp_spectra (wavelength, flux, n_points):
    lin_interp = interp1d(wavelength,flux)
    wavelength_interp = np.linspace(min(wavelength),max(flux),n_points)
    flux_interp = lin_interp(wavelength_interp)
    return wavelength_interp, flux_interp

def redden_flux (A_v, wavelength, flux):
    R_v = 3.1
    x = 1/(wavelength*10**(-4))
    y = x - 1.82

    a_x = 1+0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
    b_x = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7

    A_lambda = (a_x + (b_x/R_v))*A_v
    z = 10**(A_lambda/2.5)
    reddened_flux = flux/z

    return reddened_flux

def scale_luminosity (luminosity, flux, multiplier):
    logL = np.log10(luminosity)
    scale_1 = logL - multiplier
    scale_2 = 10.**scale_1
    scaled_flux = flux/scale_2
    return scaled_flux

def redden_scale (filename, mass, lum_multiplier, reddening):
    wavelength, flux, radius, luminosity = extract_MARCSvals (filename)
    flux_scaled = scale_luminosity(luminosity, flux, lum_multiplier)
    flux_red = redden_flux(reddening,wavelength,flux_scaled)

    return wavelength, flux_red

def combine_spectra (star1_flux, star2_flux):
    return star1_flux + star2_flux
                        

# The remaining helper functions are specific to the MARCS and BSTAR06 models
def rad_lum_from_header (filename):
    ## RSG header to determine radius
    with open(filename) as file:
        for i,line in enumerate(file):
            if i == 7:
                radius = float(line[0:12])
            if i == 8:
                luminosity = float(line[0:12])
                break
            
    return radius, luminosity

def extract_BSTAR_vals (filename):
    with open(filename,'r') as f:
        inline = f.readline().strip('\n')

    Bstar_wave, Bstar_flux = extract_BSTAR_vals(inline)
    cl = line[1]
    temperature = line[3:8]
    Mv = float(line[25:29])
    BC = float(line[33:37])
    g = line[41:44]
                
    class_dict = {"5":"V","1":"I","3":"III"}
    ty = class_dict[cl]
                
    luminosity = calc_lum(Mv,BC)
    radius = calc_radius(luminosity,float(temperature))
    
    BSTAR_model = np.loadtxt(filename,comments="#")
    wavelength, flux = interp_spectra(BSTAR_model[:,0],BSTAR_model[:,1]*radius,7000)

    return wavelength, flux

def extract_MARCS_vals (filename):
    file_flux = filename+".flx"
    file_mod = filename+".mod"

    MARCS_wave = np.loadtxt("MARCS/flx_wavelengths.vac")
    MARCS_flux = np.loadtxt(file_flux)
    MARCS_waveAir = convert_vac_to_air(MARCS_wave)

    # shorten MARCSarray
    wave_sm = MARCS_waveAir[18017:40766]
    flux_sm = MARCS_flux[18017:40766]*R_R

    # interpolate
    wavelength, flux = interp_spectra(wave_sm,flux_sm,7000)

    radius, luminosity = rad_lum_from_header(file_mod)
    
    return wavelength, flux, radius, luminosity

def main():
    # STAR 1
    # in this example, this is a red supergiant from MARCS models
    star1_filename = 's4000_g+0.0_m15._t05'
    star1_wave, star1_flux, star1_rad, star1_lum = extract_MARCS_vals (star1_filename)

    # STAR 2
    # in this example, this is a B star from the BSTAR06 models
    star2_filename = 'BG20000g400_150.11'
    star2_wave, star2_flux = extract_BSTAR_vals (star2_filename)

    # reddening range
#    reddening = np.arange(0,2.8,0.2)
    # pass in single value
    
    # mass / luminosity range
    # pass in single value
#    mass = ["10","15","20","25"]
#    lum_multiplier = [4.2,logL_R,5.2,5.4]

# np.savetxt("mods_all/RSG_"+filenameR[7:11]+"K_"+mass[i]+"Mo+B"+ty+"_"+str(T_B)+"K_R"+str(reddening[r])+".txt",np.array(list(zip(Bxnew,rFlux+Bynew))))

if __name__ == "__main__":
    main()

                        
