ACVD 
====
[![Build Status](https://travis-ci.org/valette/ACVD.png)](https://travis-ci.org/valette/ACVD)


<p align="center">
  <img src="https://www.creatis.insa-lyon.fr/~valette/public/project/acvd/featured.jpg">
</p>

### Info ###
This code is the implementation deriving from those papers:

[1] S. Valette,J.-M. Chassery and R. Prost, Generic remeshing of 3D triangular meshes with metric-dependent discrete Voronoi Diagrams, IEEE Transactions on Visualization and Computer Graphics, Volume 14, no. 2, pages 369-381, 2008.

[2] Sebastien Valette and Jean-Marc Chassery, Approximated Centroidal Voronoi Diagrams for Uniform Polygonal Mesh Coarsening, Computer Graphics Forum (Eurographics 2004 proceedings), Vol. 23, No. 3, September 2004, pp. 381-389. 

[3] M. Audette, D. Rivière, M. Ewend, A. Enquobahrie, and S. Valette, "Approach-guided controlled resolution brain meshing for FE-based interactive neurosurgery simulation", Workshop on Mesh Processing in Medical Image Analysis, in conjunction with MICCAI 2011., Toronto, Canada, pp. 176--186, 09/2011.


This code is cross-platform and should compile under Linux ,MacOS and Window$ OS.
### Licence ###
This code is distributed under the CeCILL-B license (BSD-compatible)
(copyright CNRS, INSA-Lyon, UCBL, INSERM.)


###  Dependencies ###
* VTK www.vtk.org
  * Version 6.x and 7.x : use this branch
  * Version 5.x use the vtk5 branch : https://github.com/valette/ACVD/tree/vtk5
* CMAKE www.cmake.org

###  Simple compilation HowTo under Linux ###
	git clone https://github.com/valette/ACVD.git
	cd ACVD
	cmake . -DCMAKE_BUILD_TYPE=Release
	make

the executables (ACVD, ACVDQ, AnisotropicRemeshingQ and others should be found under the "bin" subdirectory)

### Usage ###
execute ACVD and ACVDQ without arguments to see the available options.

when using graphical display, the 'e' key allows to continue to the next step during interaction

for ACVD, the output is written in the file simplification.ply

additionnally, when running ACVD, a file output_1.ply is also written. It is the output mesh before post-processing using quadrics.

note that to enforce a manifold output mesh, such as explained in [3], you need to use the -m 1 option.

comments, suggestions : https://github.com/valette/ACVD/issues
