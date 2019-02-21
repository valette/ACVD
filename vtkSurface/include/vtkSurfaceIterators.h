/*=========================================================================

  Program:   vtkSurfaceIterators
  Module:    vtkSurface
  Language:  C++
  Date:      2008/07
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

#ifndef __vtkSurfaceIterators_h
#define __vtkSurfaceIterators_h

#include "vtkSurface.h"

class vtkSurfaceVertexRingOrientedIterator
{
public:
	void SetInputData(vtkSurface *Input)
	{
		this->Input=Input;
	}

	void InitTraversal(const vtkIdType& Vertex)
	{
		this->CentralVertex=Vertex;
		vtkIdType Edge=this->Input->GetFirstEdge(Vertex);
		vtkIdType f1,f2;
		this->Input->GetEdgeFaces(Edge,f1,f2);
		vtkIdType v1,v2,v3;

		this->Input->GetEdgeVertices(Edge,this->OuterVertex,v2);

		if (this->OuterVertex==Vertex)
			this->OuterVertex=v2;

		this->Input->GetFaceVertices(f1,v1,v2,v3);
		
		if (((v1==Vertex)&&(v2==OuterVertex))
			||((v2==Vertex)&&(v3==OuterVertex))
			||((v3==Vertex)&&(v1==OuterVertex)))
			this->CurrentFace=f2;
		else
			this->CurrentFace=f1;

		this->LastVertex=this->OuterVertex;
		this->LoopFinished=false;
	}
	
	vtkIdType GetNextVertex()
	{
		if (LoopFinished)
		{
			this->CurrentFace=-1;
			this->OuterVertex=-1;
			return (-1);
		}
			
		vtkIdType NextFace,NextVertex;
		this->Input->Conquer(this->CurrentFace,this->CentralVertex,this->OuterVertex,NextFace,NextVertex);

		if ((NextVertex==this->LastVertex)&&this->OneLoopOnly)
			this->LoopFinished=true;

		this->CurrentFace=NextFace;
		this->OuterVertex=NextVertex;
		return (this->OuterVertex);
	};

	vtkIdType GetNextEdge()
	{
		return (this->Input->IsEdge(this->CentralVertex,this->GetNextVertex()));
	};
	
	vtkIdType GetVertex()
	{
		return (this->OuterVertex);
	}

	vtkIdType GetEdge()
	{
		return (this->Input->IsEdge(this->CentralVertex,this->OuterVertex));
	}

	vtkSurfaceVertexRingOrientedIterator()
	{
		this->OneLoopOnly=true;
	};

	void SetOneLoopOnly (bool Value)
	{
		this->OneLoopOnly=Value;
	}

	~vtkSurfaceVertexRingOrientedIterator(){};

private:
	vtkSurface *Input;
	bool OneLoopOnly;
	bool LoopFinished;
	vtkIdType CentralVertex;
	vtkIdType OuterVertex;
	vtkIdType CurrentFace;
	vtkIdType LastVertex;
};

class vtkSurfaceVertexRingRandomIterator
{
public:
	void SetInputData(vtkSurface *Input)
	{
		this->Input=Input;
	}

	void InitTraversal(const vtkIdType& Vertex)
	{
		this->CentralVertex=Vertex;
		this->Input->GetVertexNeighbourEdges(Vertex, this->NumberOfEdges, this->Edges);
		this->CurrentEdgeIndex=-1;
	}
	
	vtkIdType GetNextVertex()
	{
		if (++this->CurrentEdgeIndex>=this->NumberOfEdges)
			return (-1);
		this->CurrentEdge=this->Edges[this->CurrentEdgeIndex];
		return (this->GetVertex());
	};

	vtkIdType GetNextEdge()
	{
		if (++this->CurrentEdgeIndex>=this->NumberOfEdges)
			return (-1);
		this->CurrentEdge=this->Edges[this->CurrentEdgeIndex];
		return (this->CurrentEdge);
	};
	
	vtkIdType GetVertex()
	{
		vtkIdType v1,v2;
		this->Input->GetEdgeVertices(this->CurrentEdge,v1,v2);
		if (v1==this->CentralVertex)
			return (v2);
		return (v1);
	}

	vtkIdType GetEdge()
	{
		return (this->CurrentEdge);
	}

	vtkSurfaceVertexRingRandomIterator(){};
	~vtkSurfaceVertexRingRandomIterator(){};

private:
	vtkSurface *Input;
	vtkIdType CurrentEdgeIndex;
	vtkIdType CurrentEdge;
	vtkIdType *Edges;
	vtkIdType NumberOfEdges;
	vtkIdType CentralVertex;
};

#endif

