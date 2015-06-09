ACVD
====

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
* VTK www.vtk.org (version 5.x, version 6.x in the vtk6 branch : https://github.com/valette/ACVD/tree/vtk6) 
* CMAKE www.cmake.org

###  Simple compilation HowTo under Linux ###
	git clone https://github.com/valette/ACVD.git
	cd ACVD
	cmake . -DCMAKE_BUILD_TYPE=Release
	make

the executables (ACVD, ACVDQ and AnisotropicRemeshingQ should be found under the "bin" subdirectory)

### Note for window$ users ###
if you get this compilation error :
	fatal error LNK1104: cannot open file '..\bin\Debug\vtkSurface.lib'
you can fix it by unselecting 'build_shared_libs' with the cmake UI

### Options ###
execute ACVD and ACVDQ without arguments to see the available options.

the output of the program is written in the file simplification.ply

additionnally, when running ACVD, a file output_1.ply is also written. It is the output mesh before post-processing using quadrics.

note that to enforce a manifold output mesh, such as explained in [3], you need to use the -m 1 option, and to define a correct number of spare clusters using options -sf or -sc.

comments, suggestions : https://github.com/valette/ACVD/issues


