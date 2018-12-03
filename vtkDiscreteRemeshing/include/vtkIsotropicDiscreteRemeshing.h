/***************************************************************************
vtkIsotropicDiscreteRemeshing.h  -  description
-------------------
begin                : March 2006
copyright            : Sebastien Valette
email                : 
***************************************************************************/

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

#ifndef _VTKISOTROPICDISCRETEREMESHING_H_
#define _VTKISOTROPICDISCRETEREMESHING_H_

#include <vtkObjectFactory.h>

#include "vtkDiscreteRemeshing.h"
#include "vtkTrianglesProcessing.h"
#include "vtkVerticesProcessing.h"
#include "vtkIsotropicMetricForClustering.h"
#include "vtkQEMetricForClustering.h"

typedef vtkDiscreteRemeshing<vtkIsotropicMetricForClustering> IsotropicRemeshing;

/**
 * A Class to process coarsening of vtkSurface PolyData.
 * Implemented from the paper: [1] " Approximated Centroidal Voronoi Diagrams for Uniform 
 * Polygonal Mesh Coarsening", Valette & Chassery, Eurographics 2004.

 * Adaptive meshing is also possible, in spirit with:
 * [2] "Adaptive Polygonal Mesh Simplification With Discrete Centroidal Voronoi Diagrams"
 *   by, S. Valette, I. Kompatsiaris and J.-M. Chassery 
 * See the example file ACVD.CXX
 */

class VTK_EXPORT vtkIsotropicDiscreteRemeshing : public vtkVerticesProcessing<IsotropicRemeshing>
{
public:

	static vtkIsotropicDiscreteRemeshing* New()
	{
		// First try to	create the object from the vtkObjectFactory
		vtkObject* ret = vtkObjectFactory::CreateInstance("vtkIsotropicDiscreteRemeshing");
		if(ret)
		{
			return (vtkIsotropicDiscreteRemeshing*)ret;
		}
		// If the factory was unable to	create the object, then	create it here.
		return (new	vtkIsotropicDiscreteRemeshing);
	}

protected:

	vtkIsotropicDiscreteRemeshing()
	{
	}
	~vtkIsotropicDiscreteRemeshing() {};
};

typedef vtkDiscreteRemeshing<vtkQEMetricForClustering> TempQERemeshing;
class VTK_EXPORT vtkQIsotropicDiscreteRemeshing : public vtkVerticesProcessing<TempQERemeshing>
{
public:

	static vtkQIsotropicDiscreteRemeshing* New()
	{
		// First try to	create the object from the vtkObjectFactory
		vtkObject* ret = vtkObjectFactory::CreateInstance("vtkQIsotropicDiscreteRemeshing");
		if(ret)
		{
			return (vtkQIsotropicDiscreteRemeshing*)ret;
		}
		// If the factory was unable to	create the object, then	create it here.
		return (new	vtkQIsotropicDiscreteRemeshing);
	}

protected:

	vtkQIsotropicDiscreteRemeshing()
	{
	}
	~vtkQIsotropicDiscreteRemeshing() {};
};

#endif
