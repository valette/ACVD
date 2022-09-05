/*=========================================================================

  Program:   Uniform Vertices Clustering for meshes 
  Module:    vtkVerticesUniformClustering.h
  Language:  C++
  Date:      2004/09
  Author:    Sebastien VALETTE

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

#ifndef _VTKVERTICESPROCESSING_H_
#define _VTKVERTICESPROCESSING_H_

#include "vtkObjectFactory.h"

/**
 *  This class derives from vtkUniformClustering. It is aimed at constructing uniform clusters of mesh vertices
 */

template <class Processing> class  vtkVerticesProcessing : public Processing

{
public:

	static vtkVerticesProcessing *New()
	{
		vtkObject* ret = vtkObjectFactory::CreateInstance("vtkVerticesProcessing");
		if(ret)
	    {
			return (vtkVerticesProcessing*)ret;
	    }
	  // If the factory was unable to create the object, then create it here.
		return (new vtkVerticesProcessing<Processing>);
	};

	int GetNumberOfItems ()
	{
		return (this->GetInput()->GetNumberOfPoints());
	};

protected:

	vtkVerticesProcessing() 
	{
		this->ClusteringType=1;
		this->VList=vtkIdList::New();
#ifdef DOmultithread
		this->VLists=0;
#endif
	};

	~vtkVerticesProcessing()
	{
		this->VList->Delete();
#ifdef DOmultithread
		if (VLists!=0)
		{
			for (int i=0;i<this->NumberOfThreads+1;i++)
				VLists[i]->Delete();
			delete [] VQueues;
			delete [] VLists;
		}
#endif
	};

	void GetItemCoordinates(vtkIdType Item, double *P)
	{
		this->GetInput()->GetPointCoordinates(Item,P);
	}
	
	double GetItemArea(vtkIdType Item)
	{
		return (this->GetInput()->GetVertexArea(Item));
	}

	int GetNumberOfDualItems(){return (this->Input->GetNumberOfCells());};

#ifdef DOmultithread

	vtkIdList **VLists;
	std::queue<vtkIdType> *VQueues;

	void Init()
	{
		VLists=new vtkIdList*[this->NumberOfThreads+1];
		
		for (int i=0;i<this->NumberOfThreads+1;i++)
			VLists[i]=vtkIdList::New();
		VQueues=new std::queue<vtkIdType>[this->NumberOfThreads+1];
	
		this->Processing::Init();
	};

	void AddItemRingToProcess(vtkIdType Item,int ProcessId)
	{
		vtkIdType NumberOfEdges,*Edges,i;
		this->GetInput()->GetVertexNeighbourEdges(Item,NumberOfEdges,Edges);
		for (i=0;i<NumberOfEdges;i++)
			this->AddEdgeToProcess(Edges[i],ProcessId);
	};

	int ThreadedConnexityConstraintProblem(vtkIdType Item,vtkIdType Edge,vtkIdType 
Cluster, vtkIdType Cluster2, int Thread)
	{
		return ConnexityConstraintProblemLocal(Item,Edge,Cluster,this->VLists[Thread],this->VQueues[Thread]);
	};
#endif


private:

	std::queue<vtkIdType> VQueue;
	vtkIdList *VList;

	void GetItemNeighbours (vtkIdType Item, vtkIdList *IList)
	{
		this->GetInput()->GetVertexNeighbours(Item,IList); 
	};

	void AddItemRingToProcess(vtkIdType Item)
	{
		vtkIdType NumberOfEdges,*Edges,i;
		this->GetInput()->GetVertexNeighbourEdges(Item,NumberOfEdges,Edges);
		for (i=0;i<NumberOfEdges;i++)
			this->EdgeQueue.push(Edges[i]);
	};

	void GetEdgeItems(vtkIdType Item, vtkIdType &I1,vtkIdType &I2)
	{
		this->GetInput()->GetEdgeVertices(Item,I1,I2);
	};

	void GetItemEdges(vtkIdType Item, vtkIdList *EList)
	{
		this->Input->GetVertexNeighbourEdges(Item,EList);
	};

	/// check if moving the face Item from Cluster will not make Cluster non-connex
	int ConnexityConstraintProblem(vtkIdType Item,vtkIdType Edge,vtkIdType
Cluster, vtkIdType Cluster2)
	{
		return ConnexityConstraintProblemLocal(Item,Edge,Cluster,this->VList,this->VQueue);
	};

	/// check if moving the face Item from Cluster will not make Cluster non-connex
	int ConnexityConstraintProblemLocal(vtkIdType Item,vtkIdType Edge,vtkIdType
Cluster, vtkIdList *VList, std::queue<vtkIdType> &VQueue)
	{
		if (!this->ConnexityConstraint)
			return (0);
		
		VList->Reset();
		while (VQueue.size())
			VQueue.pop();

		vtkIdType v1,v2,v3,i,j,Cluster1,Cluster2;
		vtkIdType NumberOfEdges,*Edges;

		this->GetInput()->GetVertexNeighbourEdges(Item,NumberOfEdges,Edges);
		for (i=0;i<NumberOfEdges;i++)
		{
			this->Input->GetEdgeVertices(Edges[i],v1,v2);
			Cluster1=this->Clustering->GetValue(v1);
			Cluster2=this->Clustering->GetValue(v2);
			if ((v1!=Item)&&(Cluster1==Cluster))
				VList->InsertNextId(v1);
			if ((v2!=Item)&&(Cluster2==Cluster))
				VList->InsertNextId(v2);
		}

		if (VList->GetNumberOfIds()==0)
			return (0);

		vtkIdType NumberOfVertices=VList->GetNumberOfIds();
		vtkIdType NumberOfVisitedVertices=1;

		v1=VList->GetId(0);
		VList->DeleteId(v1);

		VQueue.push(v1);
		while (VQueue.size())
		{
			v1=VQueue.front();
			VQueue.pop();
			this->GetInput()->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
			for (i=0;i<NumberOfEdges;i++)
			{
                this->Input->GetEdgeVertices(Edges[i],v2,v3);
				for (j=0;j<VList->GetNumberOfIds();j++)
				{
					v1=VList->GetId(j);
                    if (v1==v2)
					{
						VQueue.push(v2);
						NumberOfVisitedVertices++;
						VList->DeleteId(v2);
						break;
					}
					else
					{
						if (v1==v3)
						{
							VQueue.push(v3);
							NumberOfVisitedVertices++;
							VList->DeleteId(v3);
							break;
						}
					}
				}
			}
		}
		if (NumberOfVertices==NumberOfVisitedVertices)
			return (0);
		return 1;
	};
};
#endif
