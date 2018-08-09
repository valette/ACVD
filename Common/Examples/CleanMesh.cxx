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

#include <vtkCommand.h>
#include "vtkSurface.h"
#include "RenderWindow.h"

int main( int argc, char *argv[] )
{
	if (argc<3)
	{
		cout<<"Usage : CleanMesh input output"<<endl;
		exit (1);
	}
	vtkSurface *Mesh=vtkSurface::New();
	Mesh->CreateFromFile(argv[1]);
	vtkPolyData *OriginalMesh=vtkPolyData::New();
	OriginalMesh->DeepCopy(Mesh);

	cout<<"Input : "<<endl;
	Mesh->DisplayMeshProperties();

	int numCells=Mesh->GetNumberOfCells();

	vtkIdList *List=vtkIdList::New();

	bool OK=false;
	int Rank=3;
	int Loop=0;
	while (!OK)
	{
		OK=true;
		int NumberOfCleanups=0;
		for (int i=0;i<numCells;i++)
		{
			if (Mesh->IsFaceActive(i))
			{
				vtkIdType vertices[4];
				Mesh->GetFaceVertices(i,vertices[0],vertices[1],vertices[2]);
				vertices[3]=vertices[0];
				vtkIdType edges[3];
				int NumberOfNonManifoldEdges=0;
				int NumberOfBoundaryEdges=0;
				for (int j=0;j<3;j++)
				{
					edges[j]=Mesh->IsEdge(vertices[j],vertices[j+1]);
					Mesh->GetEdgeFaces(edges[j],List);
					if (List->GetNumberOfIds()>2)
						NumberOfNonManifoldEdges++;

					if (List->GetNumberOfIds()==1)
						NumberOfBoundaryEdges++;
				}
				if ((NumberOfNonManifoldEdges!=0)&&(NumberOfNonManifoldEdges+NumberOfBoundaryEdges>Rank))
				{
					Mesh->DeleteFace(i);
					NumberOfCleanups++;
					OK=false;
				}
			}
		}
		cout<<"Loop "<<Loop<<" , Rank:"<<Rank<<" , "<<NumberOfCleanups<<" faces deleted"<<endl;
		Loop++;
		if (OK)
		{
			Rank--;
			OK=false;
			if (Rank==0)
				break;
		}
	}

	Mesh->CheckNormals();
	vtkSurface *CleanMesh=Mesh->CleanMemory();
	cout<<"Output : "<<endl;
	CleanMesh->DisplayMeshProperties();

	RenderWindow *Window=RenderWindow::New();
	Window->SetInputData(OriginalMesh);
	Window->Render();

	RenderWindow *Window2=RenderWindow::New();
	Window2->SetInputData(CleanMesh);
	Window2->AttachToRenderWindow(Window);
	Window2->Render();
	Window2->Interact();

	//delete objects
	Window->Delete();
	Window2->Delete();
	OriginalMesh->Delete();
	Mesh->Delete();
	CleanMesh->Delete();

	return(0);
}
