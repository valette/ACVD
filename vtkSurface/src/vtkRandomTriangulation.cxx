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

#include <vtkObjectFactory.h>
#include <vtkMath.h>
#include <vtkDelaunay2D.h>
#include "vtkRandomTriangulation.h"


//----------------------------------------------------------------------------
vtkRandomTriangulation* vtkRandomTriangulation::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkRandomTriangulation");
  if(ret) return (vtkRandomTriangulation*)ret;

  // If the factory was unable to create the object, then create it here.
  return new vtkRandomTriangulation;
}


//----------------------------------------------------------------------------
vtkRandomTriangulation::vtkRandomTriangulation()
{}

//----------------------------------------------------------------------------
vtkRandomTriangulation::~vtkRandomTriangulation()
{}

vtkSurface *vtkRandomTriangulation::BuildRandomTriangulation (int NumberOfPoints, int Type)
{
	if (Type ==7)
	{
		vtkSurface *Mesh = vtkSurface::New ();
		Mesh->AddVertex (0, 0, 0);
		Mesh->AddVertex (1, 0, 0);
		Mesh->AddVertex (0, 1, 0);
		Mesh->AddVertex (1, 1, 0);
		Mesh->AddFace (0, 1, 2);
		Mesh->AddFace (1,2,3);
		for (int i=0;i<NumberOfPoints;i++)
			Mesh->SubdivideInPlace();
		return (Mesh);
	}
	if (Type == 6)
	{
		vtkSurface *Mesh = vtkSurface::New ();
		Mesh->Init (1000, 1000, 1000);
		Mesh->AddVertex (0, 0, 0);
		Mesh->AddVertex (1, 0, 0);
		Mesh->AddVertex (0.5, 0.866, 0);
		
		Mesh->AddFace (0, 1, 2);
	
		Mesh->SubdivideInPlace();
		Mesh->SubdivideInPlace();
		Mesh->SubdivideInPlace();
		Mesh->SubdivideInPlace();
		
		Mesh->AddVertex (0.5, 0.289, 0);
		Mesh->AddVertex (0,1, 0);
		Mesh->AddVertex (1,1, 0);
		
		vtkDelaunay2D *Delaunay = vtkDelaunay2D::New ();
		vtkPolyData *Poly = vtkPolyData::New ();
		Poly->SetPoints(Mesh->GetPoints());

		Delaunay->SetInputData (Poly);
		Delaunay->Update ();
		vtkSurface *Mesh2 = vtkSurface::New ();
		Mesh2->CreateFromPolyData (Delaunay->GetOutput ());
		

		return (Mesh2);
	}
	
	if (Type == 5)
	{
		vtkSurface *Mesh = vtkSurface::New ();
		Mesh->Init (1000, 1000, 1000);
		Mesh->AddVertex (0, 0, 0);
		Mesh->AddVertex (0, 0, 1);
		Mesh->AddVertex (0, 1, 1);
		Mesh->AddVertex (0, 1, 0);
		Mesh->AddVertex (1, 0, 0);
		Mesh->AddVertex (1, 0, 1);
		Mesh->AddVertex (1, 1, 1);
		Mesh->AddVertex (1, 1, 0);

		Mesh->AddFace (0, 1, 2);
		Mesh->AddFace (0, 2, 3);
		Mesh->AddFace (1, 2, 6);
		Mesh->AddFace (1, 6, 5);
		Mesh->AddFace (0, 1, 4);
		Mesh->AddFace (1, 5, 4);
		Mesh->AddFace (4, 5, 6);
		Mesh->AddFace (4, 6, 7);
		Mesh->AddFace (0, 3, 7);
		Mesh->AddFace (0, 7, 4);
		Mesh->AddFace (3, 2, 6);
		Mesh->AddFace (3, 6, 7);

		Mesh->SubdivideInPlace ();
		Mesh->SubdivideInPlace ();
		Mesh->SubdivideInPlace ();
		Mesh->SubdivideInPlace ();
		Mesh->SubdivideInPlace ();
		Mesh->SubdivideInPlace ();

		int i;
		double P[3];
		cout << Mesh->GetNumberOfPoints () << " points" << endl;
		for (i = 0; i < Mesh->GetNumberOfPoints (); i++)
		{
			Mesh->GetPoint (i, P);
			P[0] += 0.9 * sqrt (0.6 -
					    (P[1] - 0.5) * (P[1] - 0.5));
			P[2] += 0.9 * sqrt (0.6 -
					    (P[1] - 0.5) * (P[1] - 0.5));
			Mesh->SetPointCoordinates (i, P);
		}
		return(Mesh);
	}


	int i;
	vtkPoints *Points = vtkPoints::New ();
	vtkMath *Math = vtkMath::New ();
	vtkDelaunay2D *Delaunay = vtkDelaunay2D::New ();
	vtkPolyData *Poly = vtkPolyData::New ();

	Math->RandomSeed (1000);
	float x, y, z;
	z = 0;

	if (Type < 2)
	{
		for (i = 0; i < NumberOfPoints; i++)
		{
			x = Math->Random ();
			y = Math->Random ();
			if (Type == 1)
				y = y * y * y;
			Points->InsertNextPoint (x, y, z);
		}
	}
	else
	{
		int n1, n2, n3, n4;
		n1 = 1;
		n2 = 2;
		n3 = 4;
		n4 = 8;
		Points->Allocate (NumberOfPoints * (n1 + n2 + n3 + n4));
		for (i = 0; i < n1 * NumberOfPoints; i++)
		{
			x = Math->Random ();
			y = Math->Random ();
			Points->InsertNextPoint (x, y, z);
		}
		for (i = 0; i < n2 * NumberOfPoints; i++)
		{
			x = -Math->Random ();
			y = Math->Random ();
			Points->InsertNextPoint (x, y, z);
		}
		for (i = 0; i < n3 * NumberOfPoints; i++)
		{
			x = Math->Random ();
			y = -Math->Random ();
			Points->InsertNextPoint (x, y, z);
		}
		for (i = 0; i < n4 * NumberOfPoints; i++)
		{
			x = -Math->Random ();
			y = -Math->Random ();
			Points->InsertNextPoint (x, y, z);
		}

	}

	Poly->SetPoints (Points);
	Delaunay->SetInputData (Poly);
	Delaunay->Update ();
	vtkSurface *Mesh = vtkSurface::New ();
	Mesh->CreateFromPolyData (Delaunay->GetOutput ());

	if (Type == 4)
	{
		double P[3];
		vtkSurface *Mesh2 = vtkSurface::New ();
		Mesh2->Init (10000, 10000, 10000);
		vtkIdType v1, v2, v3;

		int *OldPointsToNew = new int[Points->GetNumberOfPoints ()];
		for (i = 0; i < Points->GetNumberOfPoints (); i++)
		{
			Points->GetPoint (i, P);
			if ((P[0] < 0) && (P[1] < 0))
			{
				OldPointsToNew[i] = -1;
			}
			else
			{
				OldPointsToNew[i] = Mesh2->AddVertex (P);
			}
		}
		for (i = 0; i < Mesh->GetNumberOfCells (); i++)
		{
			Mesh->GetFaceVertices (i, v1, v2, v3);
			if ((OldPointsToNew[v1] >= 0)
			    && (OldPointsToNew[v2] >= 0)
			    && (OldPointsToNew[v3] >= 0))
			{
				Mesh2->AddFace (OldPointsToNew[v1],
						OldPointsToNew[v2],
						OldPointsToNew[v3]);
			}
		}
		Points->Delete ();
		Mesh->Delete ();
		Mesh = Mesh2;
		Points = Mesh->GetPoints ();
		//bugfix : the points are deleted in the end of the method
		Points->Register (Mesh);
	}
	if ((Type == 3) || (Type == 4))
	{
		double P[3];
		double x;
		for (i = 0; i < Points->GetNumberOfPoints (); i++)
		{
			Points->GetPoint (i, P);
			x = P[0];
			P[0] = cos (x * 1.4);
			P[2] = sin (x * 1.4);

			// here it's possible to stretch the mesh according to one direction...
			P[1] = P[1] * 4;
			Points->SetPoint (i, P);
		}
	}


	Math->Delete ();
	Points->Delete ();
	Delaunay->Delete ();
	Poly->Delete ();
	return (Mesh);
}
