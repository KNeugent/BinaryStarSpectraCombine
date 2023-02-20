import numpy as np
import matplotlib.pyplot as plt

combined = np.loadtxt("RSG_B_combined_spectrum.txt")
RSG = np.loadtxt("RSG_single_spectrum.txt")
B_star = np.loadtxt("Bstar_single_spectrum.txt")

plt.plot(RSG[:,0],RSG[:,1],'r-',label="RSG")
plt.plot(B_star[:,0],B_star[:,1],'b-',label="B star")
plt.plot(combined[:,0],combined[:,1],'k-',label="RSG+B binary")
plt.xlabel("Wavelength (Angstroms)")
plt.ylabel("Flux")
plt.legend()
plt.savefig("RSGbinarySpectrum.eps")
