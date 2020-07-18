# BioSLANT : Biodiversity on Simulated LAndscapes using Neutral Theory 

MIT landscape evolution model (Tadpole)
MIT spatially explicit neutral community model (BioSLANT)

Tadpole is a simple landscape evolution model that runs in MATLAB. The current version will run in 64-bit Mac OS X or Windows. 
BioSLANT is a spatially explicit neutral community model takes the output from Tadpole and simulates the dispersal, speciation, and extinction of organisms throughout the river networks.  

Directory: 
+ Main: Tadpole including functions for tracking and visualizing river captures and a fault with different erodibility
+ example-Tadpole fault: use one of the example scripts to get started simulating a landscape with a dipping fault of higher erodibility, calculating drainage basin reorganization statistics, and visualizing river captures.
+ BioSLANT: Neutral community model functions.
+ example-bioSLANT: Copy and paste the function "example_bioSLANt.m" into the main directory to get started. This script will walk through the main bioSLANT functions for 1) calculating the dispersal and habitat capacity functions 2) Running a NCM simulation to steady state and 3) Running BioSLANT on a dynamic landscape with river captures. Sample outputs produced by this script are in this folder as well. To use the example script, add BioSLANT and "example_Tadpole fault" to your path. 
+ StokesPerron2020_Result Example: Sample results for the fully coupled model from Stokes and Perron (in review, JGR: Earth Surface) for one set of biological parameters and file with list of parameters used in the paper. 

If you use BioSLANT for work that results in a publication, please acknowledge BioSLANT in the text and cite the following paper: 

Stokes, M.S. and Perron, J.T. (in review), Modeling the evolution of aquatic organisms in dynamic river basins. JGR: Earth Surface. 

If you use Tadpole for work that results in a publication, please acknowledge Tadpole in the text and cite one or more of the following papers, as appropriate:

Perron, J.T., W.E. Dietrich and J.W. Kirchner (2008), Controls on the spacing of ﬁrst-order valleys. Journal of Geophysical Research, 113, F04016, doi: 10.1029/2007JF000977.

Perron, J.T., J.W. Kirchner and W.E. Dietrich (2009), Formation of evenly spaced ridges and valleys. Nature, 460, 502–505, doi: 10.1038/nature08174.

Perron, J.T., P.W. Richardson, K.L. Ferrier, and M. Lapôtre (2012), The root of branching river networks. Nature, 492, 100–103, doi: 10.1038/nature11672.

Richardson, P.W., Perron, J.T., Miller, S.R., & Kirchner, J.W. (in review), Modeling the formation of topographic asymmetry by aspect-dependent erosional processes and lateral channel migration. Journal of Geophysical Research: Earth Surface. 

Richardson, P.W., Perron, J.T., Miller, S.R., & Kirchner, J.W. (in review), Unraveling the mysteries of asymmetric topography at Gabilan Mesa, California. Journal of Geophysical Research: Earth Surface.

Thanks!
