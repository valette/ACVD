
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
// Test file to show how to use the vtkUniformClustering Template Class

#include "vtkUniformClustering.h"
#include "vtkClusteringExampleFile.h"

class  vtkClusteringTest
{
public:

	vtkClusteringExample<vtkMetricExample> *Test;

			
vtkClusteringTest()
{
	Test=vtkClusteringExample<vtkMetricExample>::New();

	Test->ProcessClustering();
	Test->Delete();
}
~vtkClusteringTest()
{}

};
