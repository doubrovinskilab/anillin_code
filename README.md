Files used to simulate time-evolution of tissue invagination during ventral furrow formation in Drosophila melanogaster. All code files are thoroughly annotated throughout.

round_domain/round_domain.m
Matlab code to generate ascii files storing the initial configuration of the tissue, i.e. the coordinates of the nodes of the cellular boundaries as well as text files specifying the connectivity (i.e. which nodes are connected by edges that can carry stress). This initial configuration is the same for all simulations.
Output written to folder ./meshfiles

Other folders contain .cpp and .h files used to simulate tissue dynamics as well as matlab .m files to visualize the results. To run the code, files produced by running round_domain.m should be placed in a folder called meshfiles together with the executable. Output will be written into folder ./outpt

To visualize the output, run script cmb_.m (found in each folder with the C++ code corresponding to a given simulation). Output will be stored as .png images in folder ./mv corresponding to the simulated time-frames.
