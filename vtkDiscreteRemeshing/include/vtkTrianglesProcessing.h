/*=========================================================================

Program:   Triangles Processing (well, mostly clustering) for meshes 
Module:    vtkTrianglesProcessing.h
Language:  C++
Date:      2006/03
Auteur:    Sebastien VALETTE

=========================================================================*/

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

#ifndef _VTKTRIANGLESPROCESSING_H_
#define _VTKTRIANGLESPROCESSING_H_

#include <vtkObjectFactory.h>
#include <vtkBitArray.h>


// This class derives from vtkUniformClustering. It is aimed at constructing uniform clusters of triangles


template  <class Processing> class  vtkTrianglesProcessing : public Processing
{
public:
	
	static vtkTrianglesProcessing *New()
	{
		vtkObject* ret = vtkObjectFactory::CreateInstance("vtkTrianglesProcessing");
		if(ret)
		{
			return (vtkTrianglesProcessing*)ret;
		}
		// If the factory was unable to create the object, then create it here.
		return (new vtkTrianglesProcessing<Processing>);
	};
	
protected:


	int GetNumberOfDualItems(){return (this->Input->GetNumberOfPoints());};


	vtkTrianglesProcessing() 
	{
		this->ClusteringType=0;
		this->BoundaryVertices=0;
	}

	~vtkTrianglesProcessing()
	{
		if (this->BoundaryVertices)
			this->BoundaryVertices->Delete();
	};

	void Init()
	{
		Processing::Init();
		int i;
		if (!this->BoundaryVertices)
		{
			this->BoundaryVertices=vtkBitArray::New();
			this->BoundaryVertices->SetNumberOfValues(this->GetInput()->GetNumberOfPoints());
			for (i=0;i<this->GetInput()->GetNumberOfPoints();i++)
				if (this->GetInput()->GetNumberOfBoundaries(i)>0)
					this->BoundaryVertices->SetValue(i,1);
				else
					this->BoundaryVertices->SetValue(i,0);
		}
	}

	vtkBitArray *BoundaryVertices;

	/// returns the barycenter of the face
	void GetItemCoordinates(vtkIdType Item, double *P)
	{
		int v1,v2,v3;
		double P1[3],P2[3],P3[3];
		this->GetInput()->GetFaceVertices(Item,v1,v2,v3);
		this->GetInput()->GetPointCoordinates(v1,P1);
		this->GetInput()->GetPointCoordinates(v2,P2);
		this->GetInput()->GetPointCoordinates(v3,P3);
		P[0]=(P1[0]+P2[0]+P3[0])/3.0;
		P[1]=(P1[1]+P2[1]+P3[1])/3.0;
		P[2]=(P1[2]+P2[2]+P3[2])/3.0;
	};
	
	double GetItemArea(vtkIdType Item)
	{
		return (this->GetInput()->GetFaceArea(Item));
	}

	/// returns the list of triangles adjacent to the Item
	void GetItemNeighbours (vtkIdType Item, vtkIdList *IList)
	{
		int v1,v2,v3;
		int i;
		IList->Reset();
		vtkIdList *FList=vtkIdList::New();

		this->GetInput()->GetFaceVertices(Item,v1,v2,v3);
		this->GetInput()->GetEdgeFaces(this->GetInput()->IsEdge(v1,v2),FList);
		for (i=0;i<FList->GetNumberOfIds();i++)
			IList->InsertUniqueId(FList->GetId(i));

		this->GetInput()->GetEdgeFaces(this->GetInput()->IsEdge(v3,v2),FList);
		for (i=0;i<FList->GetNumberOfIds();i++)
			IList->InsertUniqueId(FList->GetId(i));

		this->GetInput()->GetEdgeFaces(this->GetInput()->IsEdge(v1,v3),FList);
		for (i=0;i<FList->GetNumberOfIds();i++)
			IList->InsertUniqueId(FList->GetId(i));
		FList->Delete();
	};

	/// Returns the number of faces
	int GetNumberOfItems ()
	{
		return (this->GetInput()->GetNumberOfCells());
	};

	/// adds the Item three adjacent edges to the list of edges to process
	void AddItemRingToProcess(vtkIdType Item)
	{
		int v1,v2,v3;
		this->GetInput()->GetFaceVertices(Item,v1,v2,v3);
		this->AddEdgeToProcess(this->GetInput()->IsEdge(v1,v2));
		this->AddEdgeToProcess(this->GetInput()->IsEdge(v1,v3));
		this->AddEdgeToProcess(this->GetInput()->IsEdge(v2,v3));
	};

#ifdef DOmultithread
	/// adds the Item three adjacent edges to the list of edges to process.
	/// This is the method for multithreaded clustering
	void AddItemRingToProcess(vtkIdType Item,int ProcessId)
	{
		int v1,v2,v3;
		this->GetInput()->GetFaceVertices(Item,v1,v2,v3);
		this->AddEdgeToProcess(this->GetInput()->IsEdge(v1,v2),ProcessId);
		this->AddEdgeToProcess(this->GetInput()->IsEdge(v1,v3),ProcessId);
		this->AddEdgeToProcess(this->GetInput()->IsEdge(v2,v3),ProcessId);
	};
#endif


	/// returns the two faces adjacent to the given edge
	void GetEdgeItems(vtkIdType Item, vtkIdType &I1,vtkIdType &I2)
	{
		this->GetInput()->GetEdgeFaces(Item,I1,I2);
	};

	void GetItemEdges(vtkIdType Item, vtkIdList *EList)
	{
		EList->SetNumberOfIds(3);
		int v1,v2,v3;
		this->Input->GetFaceVertices(Item,v1,v2,v3);
		EList->SetId(0,this->Input->IsEdge(v1,v2));
		EList->SetId(1,this->Input->IsEdge(v1,v3));
		EList->SetId(2,this->Input->IsEdge(v3,v2));
	}


	/// check if moving the face Item from Cluster will not make Cluster non-connex
	int ConnexityConstraintProblem(vtkIdType Item,int Edge,int Cluster)
	{
		int f1,f2,v1,v2,v3,v4,i,Cluster1;
		int NumberOfEdges,*Edges;

		if (!this->ConnexityConstraint)
			return (0);
		this->GetInput()->GetEdgeVertices(Edge,v3,v2);
		v1=this->GetInput()->GetThirdPoint(Item,v3,v2);

		this->GetInput()->Conquer(Item,v1,v2,f2,v4);
		if (f2>=0)
		{
			if (this->Clustering->GetValue(f2)!=Cluster)
				return (0);
		}
		else
			return (0);

		this->GetInput()->Conquer(Item,v1,v3,f2,v4);
		if(f2>=0)
		{
			if (this->Clustering->GetValue(f2)!=Cluster)
				return (0);
		}
		else
			return (0);

		if (this->BoundaryVertices->GetValue(v1)==1)
			return (1);

		this->GetInput()->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
		for (i=0;i<NumberOfEdges;i++)
		{
			this->GetInput()->GetEdgeFaces(Edges[i],f1,f2);
			if (f1>=0)
			{
				Cluster1=this->Clustering->GetValue(f1);
				if ((Cluster!=Cluster1)&&(Cluster1!=this->NumberOfClusters))
					return (1);
			}
			if (f2>=0)
			{
				Cluster1=this->Clustering->GetValue(f2);
				if ((Cluster!=Cluster1)&&(Cluster1!=this->NumberOfClusters))
					return (1);
			}
		}
		return (0);
	}
};
#endif
