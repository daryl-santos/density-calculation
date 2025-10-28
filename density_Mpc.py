### DENSITY CONVERSION TO MPC^-2
### Author: Daryl Joe D. Santos
### email: daryl_santos@gapp.nthu.edu.tw

from astropy.cosmology import WMAP7 # Komatsu et al. 2011 
import numpy as np
import pandas as pd

data1 = pd.read_csv('HSC_73769_subset_galaxies_density.csv')

# angular diameter distance in Mpc
z = data1['z']
d_A = WMAP7.angular_diameter_distance(z)

# get the angles in radians
theta_radian = data1['tenth_dist_corrected_subset']

# arc length = radius * angle
distance_Mpc = d_A * theta_radian
data1['tenth_dist_Mpc'] = distance_Mpc
data1['density_Mpc'] = 10/data1['tenth_dist_Mpc']**2

data1.to_csv('HSC_73769_subset_galaxies_density_Mpc.csv', index=False)
print('Done') # separation between ref. pt. and galaxy in Mpc
