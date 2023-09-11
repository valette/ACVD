/***************************************************************************
vtkNeighbourhoodComputation.h  -  description
-------------------
begin				 : January 4 2006
copyright			 : (C) 2006	by Sebastien Valette
email				 : 
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

#include <limits.h>
#include "vtkMath.h"
#include "vtkNeighbourhoodComputation.h"


void vtkNeighbourhoodComputation::ComputeNRingCells(vtkIdType Cell,int RingSize, vtkIdList *FList)
{
	this->VisitedCells->Reset();
	this->VisitedVertices->Reset();
	vtkIdType	j,k,l;
	vtkIdType v1,v2,v3,Edge1,f1,f2;
	vtkIdList *TempList;

	if (this->CellType==0)
	{
		this->Input->GetFaceVertices(Cell,v1,v2,v3);
		VList->Reset();
		VList->InsertNextId(v1);
		VList->InsertNextId(v2);
		VList->InsertNextId(v3);

		FList->Reset();
		FList->InsertNextId(Cell);
		VisitedCells->Tag(Cell);
	}
	else
	{
		VList->Reset();
		VList->InsertNextId(Cell);
		this->Input->GetVertexNeighbourFaces(Cell,FList);
		for (j=0;j<FList->GetNumberOfIds();j++)
			VisitedCells->Tag(FList->GetId(j));
	}

	vtkIdType NumberOfEdges,*Edges=0;
	for	(j=0;j<RingSize;j++)
	{
		VList2->Reset();
		for	(k=0;k<VList->GetNumberOfIds();k++)
		{	
			v1=VList->GetId(k);			
			if (!VisitedVertices->IsTagged(v1))
			{
				VisitedVertices->Tag(v1);
				this->Input->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
				for	(l=0;l<NumberOfEdges;l++)
				{
					Edge1=Edges[l];
					this->Input->GetEdgeVertices(Edge1,v2,v3);
					if (v1==v2)
						VList2->InsertNextId(v3);
					else
						VList2->InsertNextId(v2);
					this->Input->GetEdgeFaces(Edge1,f1,f2);
					if (f1>=0)
					{
						if (!VisitedCells->IsTagged(f1))
						{
							FList->InsertNextId(f1);
							VisitedCells->Tag(f1);
						}
						if (f2>=0)
						{
							if (!VisitedCells->IsTagged(f2))
							{
								FList->InsertNextId(f2);
								VisitedCells->Tag(f2);
							}
						}
					}
				}
			}
		}
		TempList=VList;
		VList=VList2;
		VList2=TempList;
	}
}
void vtkNeighbourhoodComputation::ComputeDistanceRingCells(vtkIdType Cell,double Distance, vtkIdList *FList)
{
	this->VisitedCells->Reset();
	vtkIdType v1,v2,v3,Edge1,f1,f2;
	double Point1[3],Point2[3],Point3[3],Origin[3];
	std::queue<int>	EdgesQueue;

//	if (this->CellType==0)
	{
		this->Input->GetFaceVertices(Cell,v1,v2,v3);
		VList->Reset();
		VList->InsertNextId(v1);
		VList->InsertNextId(v2);
		VList->InsertNextId(v3);

		FList->Reset();
		FList->InsertNextId(Cell);
		VisitedCells->Tag(Cell);
		EdgesQueue.push(this->GetInput()->IsEdge(v1,v2));
		EdgesQueue.push(this->GetInput()->IsEdge(v1,v3));
		EdgesQueue.push(this->GetInput()->IsEdge(v3,v2));

		this->GetInput()->GetPointCoordinates(v1,Point1);
		this->GetInput()->GetPointCoordinates(v2,Point2);
		this->GetInput()->GetPointCoordinates(v3,Point3);
		Origin[0]=(Point1[0]+Point2[0]+Point3[0])/3.0;
		Origin[1]=(Point1[1]+Point2[1]+Point3[1])/3.0;
		Origin[2]=(Point1[2]+Point2[2]+Point3[2])/3.0;

	}
	if (this->CellType!=0)
	{
		cout<<"NOT IMPLEMENTED..."<<endl;
		exit(1);
	}

	while(EdgesQueue.size())
	{
		Edge1=EdgesQueue.front();
		EdgesQueue.pop();
		this->GetInput()->GetEdgeVertices(Edge1,v1,v2);
		this->GetInput()->GetEdgeFaces(Edge1,f1,f2);
		if (VisitedCells->IsTagged(f1))
		{
			f1=f2;
		}
		if (f1>=0)
		{
			if (!VisitedCells->IsTagged(f1))
			{
				VisitedCells->Tag(f1);
				v3=this->GetInput()->GetThirdPoint(f1,v1,v2);
				this->GetInput()->GetPointCoordinates(v1,Point1);
				this->GetInput()->GetPointCoordinates(v2,Point2);
				this->GetInput()->GetPointCoordinates(v3,Point3);
				Point1[0]=(Point1[0]+Point2[0]+Point3[0])/3.0-Origin[0];
				Point1[1]=(Point1[1]+Point2[1]+Point3[1])/3.0-Origin[1];
				Point1[2]=(Point1[2]+Point2[2]+Point3[2])/3.0-Origin[2];
				FList->InsertNextId(f1);
				if (vtkMath::Norm(Point1)<Distance)
				{
					EdgesQueue.push(this->GetInput()->IsEdge(v1,v3));
					EdgesQueue.push(this->GetInput()->IsEdge(v2,v3));
				}
			}
		}
	}
}

vtkNeighbourhoodComputation::vtkNeighbourhoodComputation()
{
	this->VisitedCells=0;
	this->VisitedVertices=0;
	this->VList=vtkIdList::New();
	this->VList2=vtkIdList::New();
	this->EList=vtkIdList::New();
	this->CellType=0;
	this->VisitedCells=vtkTag::New();
	this->VisitedVertices=vtkTag::New();
}

vtkNeighbourhoodComputation::~vtkNeighbourhoodComputation()
{
	this->VList->Delete();
	this->VList2->Delete();
	this->EList->Delete();
	this->VisitedCells->Delete();
	this->VisitedVertices->Delete();
}

vtkNeighbourhoodComputation* vtkNeighbourhoodComputation::New()
{
	// First try to	create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkNeighbourhoodComputation");
	if(ret)
	{
		return (vtkNeighbourhoodComputation*)ret;
	}
	// If the factory was unable to	create the object, then	create it here.
	return (new	vtkNeighbourhoodComputation);
}


void vtkNeighbourhoodComputation::SetInputData(vtkSurface *Mesh)
{
	this->Input=Mesh;
	
	// Allocate the arrays with respect to the input mesh
	this->VisitedCells->SetNumberOfItems(Mesh->GetNumberOfCells());
	this->VisitedVertices->SetNumberOfItems(Mesh->GetNumberOfPoints());

}

