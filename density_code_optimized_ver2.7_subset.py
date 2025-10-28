### LOCAL GALAXY DENSITY CODE
### ver. 2.7
### Author: Daryl Joe D. Santos / revised by Jimmy
### email: daryl_santos@gapp.nthu.edu.tw
### Improvements:
###     1. Fixed error in edge correction
###     2. Made separate codes for normalisation and calculating density in Mpc^-2
### Major changes have been introduced. Check readme file for more info about how the code works.

import numpy as np
import pandas as pd
import math

### Define nth smallest value function
def find_nth_smallest(iterable, n, index=False):
    try:
        st = set(iterable)
        for i in range(n-1):
            st.discard(min(st))
        if index == True:
            return int(list(iterable).index(min(st)))
        else:
            return min(st)
    except: return None                     
    # Only when len(iterable) < n 

### Get data
galaxies = pd.read_csv('HSC_subset_galaxies.csv', skip_blank_lines=True)            
edgepoints = pd.read_csv('circle_edgepts.csv', skip_blank_lines=True)            # Please fill the right filename
print("The number of sources is", len(galaxies))

### Get ra, dec, and id
ra_edge = edgepoints['ra']
ra_edge = ra_edge[~np.isnan(ra_edge)]
dec_edge = edgepoints['dec']
dec_edge = dec_edge[~np.isnan(dec_edge)]
ra_gal = galaxies['ra']
ra_gal = ra_gal[~np.isnan(ra_gal)]
dec_gal = galaxies['dec']
dec_gal = dec_gal[~np.isnan(dec_gal)]
id_gal = galaxies['id']
id_gal = id_gal[~np.isnan(id_gal)]

### Convert all ra and dec to rad
ra_edge = ra_edge*(np.pi/180)
dec_edge = dec_edge*(np.pi/180)
ra_gal = ra_gal*(np.pi/180)
dec_gal = dec_gal*(np.pi/180)


### Make coordinate points/dataframe with edge and gal as x and y values
edge = list(zip(ra_edge, dec_edge))
gal = list(zip(id_gal, ra_gal, dec_gal, galaxies['z'])) # redshift values of galaxies should be named as 'redshift'
edgepd = pd.DataFrame(edge, columns=['ra','dec'])
galpd = pd.DataFrame(gal, columns=['id', 'ra','dec', 'z'])

### Append empty column to be filled later on, other columns will be created based on this column
galpd['tenth_dist_uncorrected_subset'] = ""

### Calculate sin() and cos() for each galaxy and edge first
galpd['ra_cos'] = np.cos(galpd['ra'])
galpd['dec_cos'] = np.cos(galpd['dec'])
galpd['ra_sin'] = np.sin(galpd['ra'])
galpd['dec_sin'] = np.sin(galpd['dec'])
edgepd['ra_cos'] = np.cos(edgepd['ra'])
edgepd['dec_cos'] = np.cos(edgepd['dec'])
edgepd['ra_sin'] = np.sin(edgepd['ra'])
edgepd['dec_sin'] = np.sin(edgepd['dec'])

### PART 1: Find Edge Distance
### Calculate distances from edge points for each target galaxy
def edge_calculate(galpd, edgepd):
    ### Declare empty 1D arrays where indices will be stored
    second_smallest_ind = np.empty([len(galpd)])
    smallest_ind = np.empty([len(galpd)])

    ### Declare empty arrays with ra and dec columns
    second_smallest_val = pd.DataFrame(columns=['ra','dec'])
    smallest_val = pd.DataFrame(columns=['ra','dec'])
    index = np.arange(len(galpd))

    for i, row in galpd.iterrows():
        edge_dist = np.arccos(edgepd['dec_sin']*row['dec_sin'] + edgepd['dec_cos']*row['dec_cos']*np.cos(edgepd['ra'] - row['ra']))

        ### Find indices of smallest and 2nd smallest distance
        second_smallest_ind[i] = find_nth_smallest(edge_dist, 2, index=True)
        smallest_ind[i] = find_nth_smallest(edge_dist, 1, index=True)

    ### Get nearest and 2nd nearest edgepoints as an array using the returned indices
    smallest_val['ra'] = edgepd['ra'][smallest_ind]
    smallest_val['dec'] = edgepd['dec'][smallest_ind]
    smallest_val[''] = index
    second_smallest_val['ra'] = edgepd['ra'][second_smallest_ind]
    second_smallest_val['dec'] = edgepd['dec'][second_smallest_ind]
    second_smallest_val[''] = index

    smallest_val = smallest_val.set_index('')
    second_smallest_val = second_smallest_val.set_index('')

    ### Get leg and base (angular) distance of nearest and 2nd nearest edgepoint, and galaxy
    ### Then calculate distance of galaxy from edge (altitude)
    d1 = np.arccos(np.sin(smallest_val['dec'])*galpd['dec_sin'] + np.cos(smallest_val['dec'])*galpd['dec_cos']*np.cos(smallest_val['ra'] - galpd['ra']))
    d2 = np.arccos(np.sin(second_smallest_val['dec'])*galpd['dec_sin'] + np.cos(second_smallest_val['dec'])*galpd['dec_cos']*np.cos(second_smallest_val['ra'] - galpd['ra']))
    b = np.arccos(np.sin(smallest_val['dec'])*np.sin(second_smallest_val['dec']) + np.cos(smallest_val['dec'])*np.cos(second_smallest_val['dec'])*np.cos(smallest_val['ra'] - second_smallest_val['ra']))
    s = (d1+d2+b)/2
    A = np.sqrt(s*(s-d1)*(s-d2)*(s-b))

    final_dist = (2*A/b)

    galpd['edge_dist_subset'] = final_dist # in radians
    galpd['edge_dist_arcmin_subset'] = final_dist*(60*180)/np.pi # in arcmin
    return galpd

####### AT THIS POINT, MAJOR CHANGES ARE INTRODUCED #######

### PART 2: Find 10th Nearest Neighbor Distance and Calculate Density
### Calculate distance from other galaxies for each target galaxy
### Define redshift bins for each galaxy
def distance_calculate(galpd):
    ### Declare empty arrays for better I/O
    tenth_dist_uncorrected = np.empty([len(galpd)])
    tenth_dist_val_new = np.empty([len(galpd)])
    flag = np.empty([len(galpd)])

    for i, row in galpd.iterrows():
        z = row['z']
        z_bin = (1 + z)*0.0683 # error in HSC photometry
        z_bin_upper = z + z_bin
        z_bin_lower = z - z_bin
        
        gal_arr = galpd[(galpd['z'] < z_bin_upper) & (galpd['z'] > z_bin_lower) \
                        & (galpd['ra'] != row['ra']) & (galpd['dec'] != row['dec'])]
        gal_arr = gal_arr.reset_index(drop=True)
        dist_arr = np.arccos(row['dec_sin']*gal_arr['dec_sin'] + row['dec_cos']*gal_arr['dec_cos']*np.cos(row['ra'] - gal_arr['ra']))

        ### Store distance values on temporal dataframe
        tenth_dist_uncorrected[i] = find_nth_smallest(dist_arr, 10) #in radians
        galpd.loc[i, 'tenth_dist_uncorrected_subset'] = tenth_dist_uncorrected[i]
        
### PART 3: Edge Correction
### Select galaxies near edge
        flag[i] = np.where(galpd['edge_dist_subset'][i] < tenth_dist_uncorrected[i], 1, 0)

### If the target galaxy is an edge galaxy, proceed to calculating true "10th nearest neighbor distance"
### First, identify the value of n, which corresponds to the index of the true "10th nearest neighbor"
        if flag[i] == 1:
            c = 2*np.sqrt((tenth_dist_uncorrected[i])**2 - (galpd['edge_dist_subset'][i])**2) # c = length of chord
            y = 1 - ((c**2)/(2*(tenth_dist_uncorrected[i]**2))) 
            theta = np.arccos(y) # theta = angle between 2 radii that connect on both ends of chord
            A = (tenth_dist_uncorrected[i]**2)*(theta - np.sin(theta))/2 # A = area of segment
            percent_area = A/(np.pi*(tenth_dist_uncorrected[i]**2)) # get percentage of area of segment
            n = (1 - percent_area)*10
            #n = round(n). # we round off + 1 since we do not count the target galaxy itself
            n = n.round() # For Tetsuya-san's computer
            n = int(n)
        else:
            n = 10

### Find index and then the value of 10th smallest distance
        tenth_dist_val_new[i] = find_nth_smallest(dist_arr, n)

    galpd['flag_subset'] = flag
    galpd['tenth_dist_corrected_subset'] = tenth_dist_val_new # in radians
    
    return galpd

galpd = edge_calculate(galpd, edgepd)
galpd = distance_calculate(galpd)

galpd['density_uncorrected_subset'] = 10/(galpd['tenth_dist_uncorrected_subset']**2)
galpd['tenth_dist_uncorrected_arcmin_subset'] = galpd['tenth_dist_uncorrected_subset']*(60*180)/math.pi
galpd['density_uncorrected_arcmin_subset'] = 10/(galpd['tenth_dist_uncorrected_arcmin_subset']**2)

galpd['tenth_dist_corrected_arcmin_subset'] = galpd['tenth_dist_corrected_subset']*(60*180)/math.pi
galpd['density_corrected_subset'] = 10/(galpd['tenth_dist_corrected_subset']**2)
galpd['density_corrected_arcmin_subset'] = 10/(galpd['tenth_dist_corrected_arcmin_subset']**2)

galpd.drop(['ra_cos', 'ra_sin', 'dec_cos', 'dec_sin'], axis=1, inplace=True)
galpd.to_csv('HSC_73769_subset_galaxies_density_removesort.csv', index=False)
print('Done density!')
