/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the CeCILL-B license under French law and 
*  abiding by the rules of distribution of free software. You can  use, 
*  modify and/ or redistribute the software under the terms of the CeCILL-B 
*  license as circulated by CEA, CNRS and INRIA at the following URL 
*  http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html 
*  or in the file LICENSE.txt.
*
*  As a counterpart to the access to the source code and  rights to copy,
*  modify and redistribute granted by the license, users are provided only
*  with a limited warranty  and the software's author,  the holder of the
*  economic rights,  and the successive licensors  have only  limited
*  liability. 
*
*  The fact that you are presently reading this means that you have had
*  knowledge of the CeCILL-B license and that you accept its terms.
* ------------------------------------------------------------------------ */  

// .NAME vtkRandomTriangulation - Creates meshes with random triangulation

#ifndef __VTKRANDOMTRIANGULATION
#define __VTKRANDOMTRIANGULATION

#include "vtkSurface.h"

class VTK_EXPORT vtkRandomTriangulation : public vtkObject
{
public:

	/// Description
	/// Builds a random triangulation
	/// the entry Type gives the distribution :
	/// 0 :uniform plane 
	/// 1 : non uniform plane 
	/// 2: plane made of 4 different regions with different densities
	/// 3 : half pipe
	/// 4 : half pipe cut on one corner
	static vtkRandomTriangulation *New();

	static vtkSurface * BuildRandomTriangulation (int NumberOfPoints, int Type);


protected:
  vtkRandomTriangulation();
  ~vtkRandomTriangulation();
};

#endif
