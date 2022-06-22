[![CI](https://github.com/valette/ACVD/actions/workflows/ci.yml/badge.svg)](https://github.com/valette/ACVD/actions/workflows/ci.yml)

ACVD 
====
<!---[![Build Status](https://travis-ci.org/valette/ACVD.png)](https://travis-ci.org/valette/ACVD) --> 


<p align="center">
  <img src="https://www.creatis.insa-lyon.fr/~valette/public/project/acvd/featured.jpg">
</p>

### Info ###
This code is the implementation deriving from those papers:

[1] S. Valette,J.-M. Chassery and R. Prost, Generic remeshing of 3D triangular meshes with metric-dependent discrete Voronoi Diagrams, IEEE Transactions on Visualization and Computer Graphics, Volume 14, no. 2, pages 369-381, 2008.

[2] Sebastien Valette and Jean-Marc Chassery, Approximated Centroidal Voronoi Diagrams for Uniform Polygonal Mesh Coarsening, Computer Graphics Forum (Eurographics 2004 proceedings), Vol. 23, No. 3, September 2004, pp. 381-389. 

[3] M. Audette, D. Rivière, M. Ewend, A. Enquobahrie, and S. Valette, "Approach-guided controlled resolution brain meshing for FE-based interactive neurosurgery simulation", Workshop on Mesh Processing in Medical Image Analysis, in conjunction with MICCAI 2011., Toronto, Canada, pp. 176--186, 09/2011.


This code is cross-platform and should compile under Linux, MacOS and Windows.
### Licence ###
This code is distributed under the CeCILL-B license (BSD-compatible)
(copyright CNRS, INSA-Lyon, UCBL, INSERM.)


###  Dependencies ###
* VTK www.vtk.org
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

### Multithread versions ###
For each program ACVD, ACVDQ and AnisotropicRemeshingQ, there is a parallel implementation, called ACVDP, ACVDQP and AnisotropicRemeshingQP. In the examples bellow, just add a trailing "P" to the executable to use all your processor cores. Note that the parallel versions are not deterministic, so running the programm twice with the same parameters will yield different remeshings. The parallel versions run much faster when quadrics are used (i.e. with ACVDQ or AnisotropcRemeshigQ), but the speedup is small with linear ACVD. For all programs, the number of threads can be set using the "-p numberOfThreads" option.

### Examples

#### Remeshing the Stanford bunny to 3000 vertices : ####
	wget https://github.com/alecjacobson/common-3d-test-models/raw/master/data/stanford-bunny.obj
	bin/ACVD stanford-bunny.obj 3000 0

taking into account curvature:

	bin/ACVD stanford-bunny.obj 3000 1.5

#### Remeshing the fandisk to 3000 vertices, taking into account sharp features with ACVDQ: ####
	wget https://github.com/alecjacobson/common-3d-test-models/raw/master/data/fandisk.obj
	bin/ACVDQ fandisk.obj 3000 0

#### Remeshing the horse to 1000 vertices with anisotropic metric: ####
	wget https://github.com/alecjacobson/common-3d-test-models/raw/master/data/horse.obj
	bin/AnisotropicRemeshingQ horse.obj 1000 1.5

#### Remeshing the Thai Statue to 100000 vertices with curvature computation ####
	wget http://graphics.stanford.edu/data/3Dscanrep/xyzrgb/xyzrgb_statuette.ply.gz
	gunzip xyzrgb_statuette.ply.gz
	bin/ACVDQ xyzrgb_statuette.ply 100000 1.5

	parallel version:
	bin/ACVDQP xyzrgb_statuette.ply 100000 1.5

	parallel version restricted to 3 threads:
	bin/ACVDQP xyzrgb_statuette.ply 100000 1.5 -np 3

for all the examples above, interactive visualization of the processing can be triggered by adding "-d 2" to the command lines

### Python port

A part of ACVD has been ported to python here: https://github.com/pyvista/pyacvd
