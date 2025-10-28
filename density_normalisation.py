### DENSITY NORMALISATION
### Author: Daryl Joe D. Santos
### email: daryl_santos@gapp.nthu.edu.tw
### Major changes have been introduced. Check readme file for more info about how the code works.

import numpy as np
import pandas as pd
import math

### Get data
galaxies = pd.read_csv('457183_HSC_galaxies_arcmin_corrected_radec.csv', skip_blank_lines=True) # all galaxies in NEPW field
print("The number of sources is", len(galaxies))

### Append empty columns to be filled later on
galaxies['normalised_density_corrected'] = ""
galaxies['normalised_density_corrected_arcmin'] = ""

### Define redshift bins for each galaxy
for i, row in galaxies.iterrows():
    z = row['z']
    z_bin = (1 + z)*0.0683
    z_bin_upper = z + z_bin
    z_bin_lower = z - z_bin
    gal_arr = galaxies[(galaxies['z'] < z_bin_upper) & (galaxies['z'] > z_bin_lower) & ((galaxies['ra']!=row['ra'])&(galaxies['dec']!=row['dec']))]
    gal_arr = gal_arr.reset_index(drop=True)

### Get median density per redshift bin
    median_rad = np.median(gal_arr['density_corrected'])
    median_arcmin = np.median(gal_arr['density_corrected_arcmin'])

### Normalise
    galaxies.loc[i, 'normalised_density_corrected'] = galaxies.loc[i, 'density_corrected']/median_rad
    galaxies.loc[i, 'normalised_density_corrected_arcmin'] = galaxies.loc[i, 'density_corrected_arcmin']/median_arcmin

galaxies.to_csv('457183_HSC_galaxies_arcmin_corrected_radec_normalised.csv', index=False)

print('Done!')

###### END OF CODE #######
