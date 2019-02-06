# Ab initio tomography with object heterogeneity and unknown viewing parameters
#### Arunabh Ghosh (arunabhghosh@iitb.ac.in), Ritwick Chaudhry and Ajit Rajwade

This repository contains the authors' implementation for the paper "_Ab initio tomography with object heterogeneity and unknown viewing parameters_" submitted to IEEE International Conference on Image Processing (ICIP) 2019.

 * `scripts/`: Contains our implementation of the reconstruction algorithm as described in the paper. The file `main_heterogeneity_alter.m`
 drives the entire algorithm to generate all the results as shown in the paper. All the experimental parameters such as the protein complex
 to be used, noise variance, number of conformations, etc. can be specified in the file mentioned.
  * `results/`: The reconstruction results of various sets of proteins under different levels of noise and number of conformations. 
  * `data/`: The dataset of images of different protein complexes taken from [Database of Macromolecular Movements](http://molmovdb.org/).
