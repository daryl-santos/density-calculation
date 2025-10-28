README FILE

LOCAL GALAXY DENSITY CODE 
ver. 2.7
Author: Daryl Joe D. Santos
Editor: Chih-Teng Ling
Corresponding Email: daryl_santos@gapp.nthu.edu.tw
Major changes:
	1. Faster execution time
	2. Separate codes for normalisation and converting density to Mpc^-2 for faster execution time
	3. Replaced chained assignments

Expected execution time:
	~8 hours for ~750,000 galaxies
	~30 minutes for ~80,000 galaxies

################################### EXPLANATION FOR EACH STEP ###########################################

Nth smallest value function 
	- gives the index of the nth smallest value in a 1D array
	- it discards the first n values and then takes the first element of the remaining array elements

>PART 1: Find Edge Distance
For each target galaxy, get distance from all edge points, and get the 1st and 2nd smallest distances, then use
these distances and geometry to solve for the altitude, which is the edge distance.


>PART 2: Find 10th Nearest Neighbor Distance and Calculate Density
For each target galaxy, define a subset of galaxies contained within the redshift bin defined by the target
galaxy's redshift. Then, get the distance of the target galaxy from other galaxies, then get the 10th 
smallest distance. Density is calculated by using the formula density = 10/(10th_dist)^2 where 10th_dist is
in radians (arcminutes). Unit of density is rad^-2 (arcmins^-2).

Minor change: Before, the target galaxy itself is counted, so we actually get the "11th" smallest distance
(the first smallest distance is 0, which refers to the distance of the target galaxy to itself).
Now, the target galaxy is not counted, so we change it to "10th" smallest distance.

>PART 3: Select Galaxies near Edge 
Use a flag to select galaxies near edge
If edge_dist < 10th_dist, flag = 1, galaxy is close from edge, density needs edge correction
If edge_dist > 10th_dist, flag = 0, galaxy is far from edge, density does not need edge correction

PART 4: Edge Correction
For all galaxies with flag=1, we correct their 10th nearest neighbor distance
We take into consideration the circular area outside the survey edge with 10th_dist as the radius.
If x% is the percentage area not covered by survey, we must take the nth nearest neighbor distance as the 
galaxy's true 10th nearest neighbor distance, where n = (1 - x%)*10

IMPROVEMENTS FOR FASTER EXECUTION TIME:
# 00: change range(len()) to enumerate or its equivalents
# 01: calculate sin() and cos() of every (ra, dec) at first and apeend to galpd
# 02: using numpy instead of math for array
# 03: modify find_nth_smallest() a little bit
# 04: move some array calculation out of the iteration
# 05: don't need to make a new '10th_dist_array'
# 06: optimize the vector calculation for panda.dataframe objects
# 07: new columns like 'tenth_dist_corrected' or 'edge_dist' are now added to galpd 
#     instead of galaxies dataframe for better I/O speed (you can merge them by ra, dec later)
# 08: functionalized some part for tracking performance 

#########################################################################################################

DENSITY NORMALISATION
Author: Daryl Joe D. Santos
Email: daryl_santos@gapp.nthu.edu.tw

Median density value for normalisation is not derived for each redshift bin of sources defined before 
density calculation, but is derived for the redshift bin where the source is located (same definition with 
1).

#########################################################################################################

DENSITY TO MPC^-2
Author: Daryl Joe D. Santos
Corresponding Email: daryl_santos@gapp.nthu.edu.tw

We use cosmology from Komatsu et al. 2011:
	H_0 = 70.4 km/s/Mpc
	Om_L = 0.728
	Om_M = 0.272
tenth_dist_corrected is converted to Mpc first: I = tenth_dist_corrected [rad] * angular_distance [Mpc, depends on redshift of source]

#########################################################################################################

# note: why don't parallelize them? 
'''def parallelize_dataframe(df, func, n_cores=4):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df'''


