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

#include <stack>
#include <assert.h>
#include <vtkObjectFactory.h>
#include <vtkMath.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkIdListCollection.h>

#include "vtkSurfaceBase.h"

void vtkSurfaceBase::DisplayInternals( bool exitOnError ) {

	cout << "********************************************" << endl;
	cout << "*******    Internals:                  *****" << endl;
	cout << this->GetNumberOfPoints() << " Vertices " << endl;
	cout << this->GetNumberOfCells() << " Polygons " << endl;
	cout << this->GetNumberOfEdges() << " Edges " << endl;
	cout << "********************************************" << endl;
	cout << "*******    Polygons:                  ******" << endl;
	for ( int i = 0; i < this->GetNumberOfCells(); i++ ) {
		vtkIdType nVertices, *vertices;
		cout << "Face " << i << " vertices :";
		this->GetFaceVertices( i, nVertices, vertices );
		for ( int j = 0; j < nVertices; j++ ) cout << " " << vertices[ j ];
		cout << endl;
	}
	vtkIdList *list = vtkIdList::New();
	cout << "********************************************" << endl;
	cout << "*******    Edges:                     ******" << endl;
	for ( int i = 0; i < this->GetNumberOfEdges(); i++ ) {
		Edge &e = this->Edges[ i ];
		cout << "Edge " << i << " vertices :";
		cout << e.Vertex1 << " " << e.Vertex2;
		this->GetEdgeFaces( i, list );
		if ( list->GetNumberOfIds() ) {
			cout << ", faces :";
			for ( int j = 0; j < list->GetNumberOfIds(); j++)
				cout << " " << list->GetId( j );
		} else cout << ", no adjacent face";
		cout << endl;
	}

	list->Delete();
	bool error = this->CheckStructure();
	if ( error && exitOnError ) exit( 1 );
}

int vtkSurfaceBase::GetEdgeNumberOfAdjacentFaces(const vtkIdType &e)
{
	vtkIdList *NonManifoldFaces=this->Edges[e].NonManifoldFaces;
	if (NonManifoldFaces!=0)
		return (2+NonManifoldFaces->GetNumberOfIds());

	vtkIdType f1,f2;
	this->GetEdgeFaces(e,f1,f2);
	if (f2!=-1)
		return (2);
	if (f1!=-1)
		return (1);
	else
		return (0);
}

vtkIdType vtkSurfaceBase::BisectEdge(vtkIdType e)
{
	vtkIdType v1,v2,v3;
	vtkIdList *FList=vtkIdList::New();
	this->GetEdgeVertices(e,v1,v2);
	this->GetEdgeFaces(e,FList);
	
	double P1[3];
	double P2[3];

	this->GetPoint(v1,P1);
	this->GetPoint(v2,P2);
	for (int i=0;i<3;i++)
		P1[i]=0.5*(P1[i]+P2[i]);

	vtkIdType NewVertex=this->AddVertex(P1);
	this->DeleteEdgeInRing(e,v1);
	Edge &edge = this->Edges[ e ];
	edge.Vertex1 = NewVertex;
	edge.Vertex2 = v2;

	for (int i=0;i<FList->GetNumberOfIds();i++)
	{
		vtkIdType Face=FList->GetId(i);
		v3=this->GetThirdPoint(Face,v1,v2);
		this->DeleteFaceInRing(Face,this->IsEdge(v1,v3));
		
		vtkIdType ve1,ve2,ve3;
		this->GetFaceVertices(Face,ve1,ve2,ve3);
		if (((ve1==v1)&&(ve2==v2))
			||((ve2==v1)&&(ve3==v2))
			||((ve3==v1)&&(ve1==v2)))
		{
			this->AddFace(v1,NewVertex,v3);
			this->SetFace(Face,NewVertex,v2,v3);
		}
		else
		{
			this->AddFace(v3,NewVertex,v1);
			this->SetFace(Face,NewVertex,v3,v2);
		}

		this->AddEdge(v3,NewVertex,Face);
	}
	
	this->InsertEdgeInRing(e,NewVertex);

	FList->Delete();
	return (NewVertex);
}

vtkIdType FindVertexIndex (const vtkIdType *Vertices,const vtkIdType &Vertex, const vtkIdType &NumberOfVertices)
{
	for (vtkIdType i=0;i<NumberOfVertices;i++)
	{
		if (Vertex==Vertices[i])
			return (i);
	}
	return (-1);
}

void vtkSurfaceBase::ChangeFaceVertex(vtkIdType Face, vtkIdType OldVertex, vtkIdType NewVertex)
{
	vtkIdType *Vertices;
	vtkIdType NumberOfVertices;
	this->GetFaceVertices(Face,NumberOfVertices,Vertices);

	vtkIdType Index=FindVertexIndex(Vertices,OldVertex,NumberOfVertices);
	assert (Index>=0);

	vtkIdType Neighbours[2];
	Neighbours[0]=Vertices[(Index+1)%NumberOfVertices];
	Neighbours[1]=Vertices[(Index+NumberOfVertices-1)%NumberOfVertices];

	Vertices[Index]=NewVertex;

	for (int i=0;i<2;i++)
	{
		vtkIdType e=this->IsEdge(OldVertex,Neighbours[i]);
		Edge &edge = this->Edges[ e ];
		if ((edge.Poly2<0)&&(this->IsEdge(NewVertex,Neighbours[i])<0))
		{
			// this edge has only one adjacent face (the face we are modifying right now).
			// we just modify the edge
			edge.Vertex1 = NewVertex;
			edge.Vertex2 = Neighbours[i];
			this->DeleteEdgeInRing(e,OldVertex);
			this->InsertEdgeInRing(e,NewVertex);
			this->CleanVertex(OldVertex);
		}
		else
		{
			// this edge has other adjacent faces, we will create a new one.
			this->DeleteFaceInRing(Face,e);
			this->CleanEdge(e);
			this->AddEdge(NewVertex,Neighbours[i],Face);
		}
	}

	VisitedPolygons->SetValue(Face,0);
	this->ConquerOrientationFromFace(Face);
	this->Polys->Modified();
}

void vtkSurfaceBase::SetOrientationOn()
{
	this->OrientedSurface=true;
}

void vtkSurfaceBase::SetOrientationOff()
{
	this->OrientedSurface=false;
}

bool vtkSurfaceBase::CheckStructure()
{
	bool Problem=false;
	vtkIdType *Vertices,NumberOfVertices;

	for (vtkIdType i=0;i<this->GetNumberOfCells();i++)
	{
		if (this->IsFaceActive(i))
		{
			this->GetFaceVertices(i,NumberOfVertices,Vertices);
			for (vtkIdType j=0;j<NumberOfVertices;j++)
			{
				if (this->IsEdge(Vertices[j],Vertices[(j+1)%NumberOfVertices])<0)
				{
					cout<<"Problem ! Face "<<i<<" misses edge ["<<Vertices[j]<<" "
					<<Vertices[(j+1)%NumberOfVertices]<<"]"<<endl;
					Problem=true;
				}
			}
		}
	}

	vtkIdList *fList = vtkIdList::New();
	for ( int e = 0; e < this->Edges.size(); e++) {

		this->GetEdgeFaces( e, fList );
		vtkIdType v1, v2;
		this->GetEdgeVertices( e, v1, v2 );
		vtkIdType v[] = { v1, v2 };

		for ( int i = 0; i < fList->GetNumberOfIds(); i++ ) {

			vtkIdType Face = fList->GetId( i );
			vtkIdType nVertices, *vertices;
			this->GetFaceVertices( Face, nVertices, vertices );
			for ( int j = 0; j < 2; j++ ) {
				bool found = false;
				for ( int k = 0; k < nVertices; k++ )
					if ( v[ j ] == vertices[ k ] ) found = true;

				if ( !found ) {
					cout << "Problem ! In edge " << e <<": Vertex "
						<< v[ j ] << " should be in face " << Face
						<< " but it is not." << endl;
					Problem=true;
				}
			}
		}

	}

	fList->Delete();
	return Problem;
}

bool vtkSurfaceBase::IsVertexManifold( const vtkIdType& iV )
{
	vtkIdType v1,v2;
	vtkIdType FirstEdge,FirstVertex;
	vtkIdType f1,f2,f3;
	vtkIdType *Edges,NumberOfRemainingEdges;
	
	this->GetVertexNeighbourEdges(iV,NumberOfRemainingEdges,Edges);

	// if there is only one edge, the vertex is non manifold
	if (NumberOfRemainingEdges<2)
		return (false);
		
	// detect non-manifold edges
	for (vtkIdType i=0;i<NumberOfRemainingEdges;i++)
	{
		if (!this->IsEdgeManifold(Edges[i]))
			return (false);
	}

	FirstEdge=this->GetFirstEdge(iV);
	this->GetEdgeVertices(FirstEdge,FirstVertex,v1);
	if (FirstVertex==iV)
		FirstVertex=v1;
	v1=FirstVertex;

	NumberOfRemainingEdges--;
	this->GetEdgeFaces(FirstEdge,f1,f2);
	
	// turn in the first direction
	v2=this->GetThirdPoint(f1,iV,FirstVertex);
	do
	{
		if (--NumberOfRemainingEdges==0)
			return (true);
		Conquer(f1,iV,v2,f3,v1);
		f1=f3;
		v2=v1;
	} while ((f1>=0)&&(v2!=FirstVertex));

	// if we turned around back to the first vertex, or there is only one direction (f2==-1) return false
	if ((f2<0)||(v2==FirstVertex))
		return (false);

	// turn in the second direction
	v1=FirstVertex;
	v2=this->GetThirdPoint(f2,iV,FirstVertex);
	do
	{
		if (--NumberOfRemainingEdges==0)
			return (true);
		Conquer(f2,iV,v2,f3,v1);
		f2=f3;
		v2=v1;
		
	} while ((f2>=0)&&(v2!=FirstVertex));
	// all adjacent faces have been visited, but there are edges remaining, so return false
	return (false);
}

/// switches the cells orientation (usefull when the mesh is displayed all black...)
void vtkSurfaceBase::SwitchOrientation()
{
	vtkIdType NumberOfVertices, *Vertices;
	this->GetFaceVertices(0,NumberOfVertices,Vertices);
	vtkIdType *Vertices2=new vtkIdType[NumberOfVertices];
	vtkIdType i;
	for (i=0;i<NumberOfVertices;i++)
		Vertices2[i]=Vertices[i];

	for (i=0;i<NumberOfVertices;i++)
		Vertices[i]=Vertices2[NumberOfVertices-1-i];

	this->CheckNormals();	
	delete [] Vertices2;
}

void vtkSurfaceBase::SQueeze()
{

	// Resize vertices Atributes
	vtkIdType numVertices=this->GetNumberOfPoints();
	this->VerticesAttributes.resize(numVertices);
	this->ActiveVertices->Resize(numVertices);

	// Resize polygons Atributes
	this->VisitedPolygons->Resize(this->GetNumberOfCells());

	// Squeeze the PolyData
	this->vtkPolyData::Squeeze();
	this->Polys->Modified();
}

void vtkSurfaceBase::ConquerOrientationFromFace(vtkIdType Face)
{
	std::queue <vtkIdType> EdgesQueue;
	vtkIdType v1,v2,v3,v4,f1,f2,Edge1,Edge2;
	vtkIdType Visited1,Visited2;
	vtkIdType NumberOfPoints1,NumberOfPoints2,*Face1,*Face2;
	vtkIdType j;
	vtkIdType FlipFace;

	this->GetFaceVertices(Face,NumberOfPoints1,Face1);
	for (j=0;j<NumberOfPoints1;j++)
	{
		Edge1=this->IsEdge(Face1[j],Face1[(j+1)%NumberOfPoints1]);
		if (this->IsEdgeManifold(Edge1))
			EdgesQueue.push(Edge1);
	}

	int Index1,Index2;

	while (EdgesQueue.size())
	{
		// pop an edge from the queue
		Edge1=EdgesQueue.front();
		EdgesQueue.pop();
		this->GetEdgeFaces(Edge1,f1,f2);
		if (f2>=0)
		{
			Visited1=this->VisitedPolygons->GetValue(f1);
			Visited2=this->VisitedPolygons->GetValue(f2);
			if ((Visited1!=Visited2))
			{
				// if one adjacent face was not visited yet
				if (Visited2==1)
				{
					v1=f1;
					f1=f2;
					f2=v1;
					Visited1=Visited2;
				}
				// f1 was wisited, not f2
				FlipFace=0;
				this->GetFaceVertices(f1,NumberOfPoints1,Face1);
				this->GetFaceVertices(f2,NumberOfPoints2,Face2);
				this->GetEdgeVertices(Edge1,v1,v2);
				Index1=FindVertexIndex(Face1,v1,NumberOfPoints1);
				Index2=FindVertexIndex(Face2,v1,NumberOfPoints2);

				v3=Face1[(Index1+1)%NumberOfPoints1];
				v4=Face2[(Index2+1)%NumberOfPoints2];
				if (((v3==v2)&&(v4==v2))||((v3!=v2)&&(v4!=v2)))
				{
					// flip f2
					for (j=0;j<NumberOfPoints2/2;j++)
					{
						v1=Face2[j];
						Face2[j]=Face2[NumberOfPoints2-1-j];
						Face2[NumberOfPoints2-1-j]=v1;
					}
				}

				this->VisitedPolygons->SetValue(f2,1);
				// add adjacent edges of the face to the queue, except Edge
				for (j=0;j<NumberOfPoints2;j++)
				{
					Edge2=this->IsEdge(Face2[j],Face2[(j+1)%NumberOfPoints2]);
					if ((Edge2!=Edge1)&&(this->IsEdgeManifold(Edge2)==1))
						EdgesQueue.push(Edge2);
				}
			}
		}
	}
}

double vtkSurfaceBase::GetValenceEntropy()
{

	int i;		// Loop counter
	int v1;	// Valence of given point
	double inv_number;	// 1.0 / number of points
	double inv_Log2;	// 1.0 / log(2.0)
	double s;			// 
	double p;			//

	// Array of valences, indexed by point's ID
	vtkIntArray *Vals = vtkIntArray::New();
	Vals->Resize(1000);

	// Compute valence for each point and fill Vals
	vtkIdList *List = vtkIdList::New();
	for (i=0; i<1000; i++) Vals->SetValue(i,0);
	for (i=0; i<this->GetNumberOfPoints(); i++)
	{
		this->GetVertexNeighbourFaces(i,List);
		v1 = List->GetNumberOfIds();
		Vals->SetValue(v1,Vals->GetValue(v1)+1);

	}
	List->Delete();

	inv_number = 1.0 / this->GetNumberOfPoints();
	inv_Log2 = 1.0 / log(2.0);
	s = 0;

	for (i=0;i<1000;i++)
	{
		p = Vals->GetValue(i);
		if (p)
		{
			p *= inv_number;
			s -= p * log(p) * inv_Log2;
		}
	}

	Vals->Delete();

	return (s);
}

void vtkSurfaceBase::GetValenceTab(char *filename)
{
	int i;
	int v1;

	int min=0;
	int Max=0;

	// allocate array for valences
	vtkIntArray *Vals=vtkIntArray::New();
	Vals->Resize(1000);
	for (i=0;i<1000;i++) Vals->SetValue(i,0);

	// fill valence's array
	vtkIdList *List=vtkIdList::New();
	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->GetVertexNeighbourFaces(i,List);
		v1=List->GetNumberOfIds();
		Vals->SetValue(v1,Vals->GetValue(v1)+1);
	}
	List->Delete();


	// find first non zero value's index
	for(i=0;i<1000;i++)
	{
		if((Vals->GetValue(i)==0)&&(Vals->GetValue(i+1)!=0))
		{
			min = i;
			break;
		}
	}

	// find last non zero value's index
	for(i=0;i<1000;i++)
	{
		if((Vals->GetValue(i)!=0)&&(Vals->GetValue(i+1)==0))
			Max = i+2;
	}

	cout <<min <<endl;
	cout <<Max <<endl;

	std::ofstream out;
	out.open (filename, std::ofstream::out | std::ofstream::trunc);

	for(i=min;i<Max;i++)
	{
		if(Vals->GetValue(i)!=0)
			out<<"Degree "<<i<<" : "<<Vals->GetValue(i)+1<<endl;
		else
			out<<"Degree "<<i<<" : 0"<<endl;
	}

	out.close();

	Vals->Delete();

}

void vtkSurfaceBase::GetEdgeFaces(vtkIdType e1,vtkIdList *FList)
{
	FList->Reset();
	Edge &edge = this->Edges[ e1 ];
	vtkIdType Face=edge.Poly1;
	if (Face==-1)
		return;
	FList->InsertNextId(Face);
	Face=edge.Poly2;
	if (Face==-1)
		return;
	FList->InsertNextId(Face);
	vtkIdList *FList2=edge.NonManifoldFaces;
	if (FList2==0)
		return;
	int NumberOfFaces=FList2->GetNumberOfIds()-1;
	for (;NumberOfFaces!=-1;NumberOfFaces--)
		FList->InsertNextId(FList2->GetId(NumberOfFaces));
}

/****************************************************************/
/****************************************************************/
void vtkSurfaceBase::DeleteVertex(vtkIdType v1)
{
	this->VerticesGarbage.push(v1);
	this->ActiveVertices->SetValue(v1,0);
}
// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::DeleteEdge(vtkIdType EdgeToRemove)
{
	Edge &e = this->Edges[ EdgeToRemove ];
	if (e.Poly1!=-1)
		cout<<"ERROR: Trying to remove an edge which is not free !"<<endl;

	e.Active=false;
	this->EdgesGarbage.push(EdgeToRemove);
	vtkIdType v1,v2;
	this->GetEdgeVertices(EdgeToRemove,v1,v2);
	this->DeleteEdgeInRing(EdgeToRemove,v1);
	this->DeleteEdgeInRing(EdgeToRemove,v2);
	this->CleanVertex(v1);
	this->CleanVertex(v2);
	e.Vertex2 = v1;
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::DeleteFace(vtkIdType f1)
{
	// test whether the face was already deleted
	if (this->IsFaceActive(f1)==0) return;
		
	vtkIdType NumberOfPoints,*Points;
	this->GetFaceVertices(f1,NumberOfPoints,Points);

	for (vtkIdType i=0;i<NumberOfPoints;i++) {
		vtkIdType e=this->IsEdge(Points[i],Points[(i+1)%NumberOfPoints]);
		this->DeleteFaceInRing(f1,e);
		this->CleanEdge(e);
	}

	for (vtkIdType i=0;i<NumberOfPoints;i++)
		Points[i]=Points[0];

	#if ( (VTK_MAJOR_VERSION < 9))
//	this->DeleteCell(f1);
	#endif
	this->CellsGarbage[NumberOfPoints].push(f1);
	this->ActivePolygons->SetValue(f1,0);
	this->Polys->Modified();
	this->Modified();
}
void vtkSurfaceBase::CleanEdge(const vtkIdType &e)
{
	if (this->CleanEdges==0) return;
	if (this->Edges[e].Poly1==-1) this->DeleteEdge(e);
}
void vtkSurfaceBase::CleanVertex(const vtkIdType &Vertex)
{
	if (this->CleanVertices==0)
		return;
	if (this->GetValence(Vertex)==0)
		this->DeleteVertex(Vertex);
}


void vtkSurfaceBase::MergeVertices(vtkIdType v1, vtkIdType v2)
{
	vtkIdList *List=vtkIdList::New();
	vtkIdList *List2=vtkIdList::New();
	
	vtkIdType e=this->IsEdge(v1,v2);
	
	// delete polygons adjacent to the edge [v1 v2] (if there are any)
	if (e>=0)
	{
		this->GetEdgeFaces(e,List);
		for (int i=0;i<List->GetNumberOfIds();i++)
			this->DeleteFace(List->GetId(i));
	}
	
	// for each remaining face adjacent to v2, replace v2 by v1
	this->GetVertexNeighbourFaces(v2,List);
	for (vtkIdType i=0;i<List->GetNumberOfIds();i++)
	{
		vtkIdType NumberOfVertices;
		vtkIdType *Vertices;
		GetFaceVertices(List->GetId(i), NumberOfVertices, Vertices);
		Vertices[FindVertexIndex (Vertices,v2,NumberOfVertices)]=v1;
	}
	
	// Modify every edge adjacent to v2
	this->GetVertexNeighbourEdges(v2,List);
	for (vtkIdType i=0;i<List->GetNumberOfIds();i++)
	{
		vtkIdType v3,v4;
		e=List->GetId(i);
		this->GetEdgeVertices(e,v3,v4);
		if (v4==v2)
		{
			// ensure that v3==v2
			vtkIdType Vertex=v3;
			v3=v4;
			v4=Vertex;
		}
		
		vtkIdType Edge2;
		Edge2=this->IsEdge(v1,v4);
		if (Edge2>=0)
		{
			// the edge [v1 v4] already exists. Merge it with [v2 v4]
			this->GetEdgeFaces(e,List2);
			for (int j=0;j<List2->GetNumberOfIds();j++)
			{
				vtkIdType Face;
				Face=List2->GetId(j);
				this->DeleteFaceInRing(Face,e);
				this->InsertFaceInRing(Face,Edge2);
			}
			
			// delete [v2 v4]
			this->DeleteEdge(e);
		}
		else
		{
			// the edge [v1 v4] does not exist. Just modify [v2 v4] to [v1 v4]
			Edge &edge = this->Edges[ e ];
			edge.Vertex1 = v1;
			edge.Vertex2 = v4;

			this->DeleteEdgeInRing(e,v2);
			this->InsertEdgeInRing(e,v1);
		}
	}
	
	this->DeleteVertex(v2);
	List->Delete();
	List2->Delete();
	this->Polys->Modified();
}

vtkSurfaceBase* vtkSurfaceBase::New()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkSurfaceBase");
	if(ret)
	{
		return (vtkSurfaceBase*)ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkSurfaceBase);
}

void vtkSurfaceBase::CheckNormals()
{
	int i;
	for (i=0;i<this->GetNumberOfCells();i++)
		this->VisitedPolygons->SetValue(i,0);

	// We have to loop on all the mesh faces to handle correctly non connnected components
	for (i=0;i<this->GetNumberOfCells();i++)
	{
		if ((this->VisitedPolygons->GetValue(i)==0)&&(this->IsFaceActive(i)==1))
		{
			this->VisitedPolygons->SetValue(i,1);
			this->ConquerOrientationFromFace(i);
		}
	}
	this->OrientedSurface=true;
}

vtkIdType vtkSurfaceBase::IsFace(const vtkIdType &v1, const vtkIdType &v2, const vtkIdType &v3)
{
	vtkIdType edge;
	vtkIdType f1;
	vtkIdType f2;
	int F;
	edge = this->IsEdge(v1,v2);
	
	// test if the edge [v1v2] exists or not
	if (edge<0) 
		return (-1);

	this->GetEdgeFaces(edge,f1,f2);
	
	// test whether the edge [v1v2] is isolated or not
	if (f1<0) 
		return (-1);
		
	if (v3 == this->GetThirdPoint(f1,v1,v2)) 
		return (f1);
	
	// test whether if there is another adjacent face
	if (f2<0) 
		return (-1);

	if (v3==this->GetThirdPoint(f2,v1,v2)) 
		return (f2);

	// test whether there are non-manifold adjacent faces
	vtkIdList *FList2=this->Edges[edge].NonManifoldFaces;
	if (FList2==0)
		return (-1);
	
	int NumberOfFaces=FList2->GetNumberOfIds()-1;
	for (;NumberOfFaces!=-1;NumberOfFaces--)
	{
		F=FList2->GetId(NumberOfFaces);
		if (v3==this->GetThirdPoint(F,v1,v2))
		return (F);
    }
	return (-1);
}

void vtkSurfaceBase :: GetFaceNeighbours(vtkIdType Face,vtkIdListCollection *FList)
{
	vtkIdType v0,v1,v2;
	vtkIdType e0,e1,e2;
		
	vtkIdList *e0Face = vtkIdList :: New();
	vtkIdList *e1Face = vtkIdList :: New();
	vtkIdList *e2Face = vtkIdList :: New();
	
	this->GetFaceVertices(Face,v0,v1,v2);
	e0=this->IsEdge(v0,v1);
	e1=this->IsEdge(v0,v2);
	e2=this->IsEdge(v1,v2);
		
	if(e0!=-1)
	{
		this->GetEdgeFaces(e0,e0Face);
		FList->AddItem(e0Face);
	}

	if(e1!=-1)
	{
		this->GetEdgeFaces(e1,e1Face);
		FList->AddItem(e1Face);
	}

	if(e2!=-1)
	{
		this->GetEdgeFaces(e2,e2Face);
		FList->AddItem(e2Face);
	}
		
	e0Face->Delete();
	e1Face->Delete();
	e2Face->Delete();
}

// ****************************************************************
// ****************************************************************

void vtkSurfaceBase::GetFaceNeighbours(vtkIdType Face,vtkIdList *FList)
{
	vtkIdType NumberOfVertices, *Vertices;
	vtkIdList *List = vtkIdList :: New();
	this->GetFaceVertices(Face, NumberOfVertices, Vertices);
	FList->Reset();

	for (vtkIdType  i=0;i<NumberOfVertices;i++)
	{
		vtkIdType e;
		if (i<NumberOfVertices-1)
			e=this->IsEdge(Vertices[i],Vertices[i+1]);
		else
			e=this->IsEdge(Vertices[i],Vertices[0]);

		this->GetEdgeFaces(e,List);
		for (vtkIdType j=0;j<List->GetNumberOfIds();j++)
			FList->InsertUniqueId(List->GetId(j));
	}
	List->Delete();
}

vtkIdType vtkSurfaceBase::FlipEdgeSure(vtkIdType edge)
{
	std::stack<vtkIdType> Edges;
	Edges.push(edge);
	while (!Edges.empty())
	{
		vtkIdType EdgeToFlip=Edges.top();
		vtkIdType ExistingEdge=this->FlipEdge(EdgeToFlip);
		if (ExistingEdge!=-1)
			Edges.push(ExistingEdge);
		else
			Edges.pop();
	}
	return (-1);
}

// ****************************************************************
// ****************************************************************
 vtkIdType vtkSurfaceBase::FlipEdge(vtkIdType edge)
{
	vtkIdType edge1,v1,v2,v3,v4,v7,v8,v9,f1,f2;

	this->GetEdgeFaces(edge,f1,f2);
	if (f2<0)
	{
//		cout<<"*** Problem: flip edge <<"<<edge<<" which is not adjacent to 2 faces !***"<<endl;
		return(edge);
	}
	this->GetEdgeVertices(edge,v1,v2);

	v3=this->GetThirdPoint(f1,v1,v2);
	v4=this->GetThirdPoint(f2,v1,v2);
	edge1=this->IsEdge(v3,v4);
	if (edge1>=0)
	{
//		cout<<"*** Problem: flip edge "<<edge<<" but the resulting edge already exists!***"<<endl;
		return (edge1);
	}
	this->DeleteEdgeInRing(edge,v1);
	this->DeleteEdgeInRing(edge,v2);

	Edge &e = this->Edges[ edge ];
	e.Vertex1 = v3;
	e.Vertex2 = v4;

	this->InsertEdgeInRing(edge,v3);
	this->InsertEdgeInRing(edge,v4);

	this->GetFaceVertices(f1,v7,v8,v9);
	if (v1==v7)
	{
		if (v2==v8)
		{
			this->SetFace(f1,v7,v4,v9);
			this->SetFace(f2,v4,v8,v9);
		}
		else
		{
			this->SetFace(f1,v7,v8,v4);
			this->SetFace(f2,v4,v8,v9);
		}
	}
	else
	{
		if (v1==v8)
		{
			if (v2==v9)
			{
				this->SetFace(f1,v8,v4,v7);
				this->SetFace(f2,v7,v4,v9);
			}
			else
			{
				this->SetFace(f1,v9,v4,v8);
				this->SetFace(f2,v4,v9,v7);
			}
		}
		else
		{
			if (v2==v7)
			{
				this->SetFace(f1,v1,v4,v8);
				this->SetFace(f2,v8,v4,v7);
			}
			else
			{
				this->SetFace(f1,v9,v3,v4);
				this->SetFace(f2,v8,v4,v3);
			}
		}
	}

	edge1=this->IsEdge(v1,v4);
	Edge &e1 = this->Edges[ edge1 ];
	if (e1.Poly1==f2)
		e1.Poly1 = f1;
	else
		e1.Poly2 = f1;

	edge1=this->IsEdge(v2,v3);
	Edge &e2 = this->Edges[ edge1 ];
	if (e2.Poly1==f1)
		e2.Poly1 =f2;
	else
		e2.Poly2= f2;

	this->Polys->Modified();
	return (-1);
}

// ** METHODE GetNumberOfBoundaries
int vtkSurfaceBase::GetNumberOfBoundaries(const vtkIdType &v1)
{
	int numBoundaries = 0;
	vtkIdType *Edges,NumberOfEdges;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);

	for (vtkIdType i=0;i<NumberOfEdges;i++)
		if (this->Edges[i].Poly2<0)	numBoundaries++;

	return (numBoundaries/2);
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::GetVertexNeighbourEdges(vtkIdType v1, vtkIdList *Output)
{
	Output->Reset();
	vtkIdType *Edges,NumberOfEdges;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
	for (vtkIdType i=0;i<NumberOfEdges;i++)
		Output->InsertNextId(Edges[i]);
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::GetVertexNeighbourFaces(const vtkIdType &v1, vtkIdList *Output)
{
	vtkIdType f1;
	vtkIdType f2;
	vtkIdType NumberOfEdges,*Edges;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
	Output->Reset();
	for (int i=0;i<NumberOfEdges;i++)
	{
		this->GetEdgeFaces(Edges[i],f1,f2);
		if (f1>=0) Output->InsertUniqueId(f1);
		if (f2>=0) Output->InsertUniqueId(f2);
		
		vtkIdList *OtherFaces=this->Edges[Edges[i]].NonManifoldFaces;
		if (OtherFaces)
		{
			for (int j=0;j<OtherFaces->GetNumberOfIds();j++)
				Output->InsertUniqueId(OtherFaces->GetId(j));
		}
	}
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::GetVertexNeighbours(vtkIdType v1, vtkIdList *Output)
{
	Output->Reset();
	vtkIdType NumberOfEdges,*NeighbourEdges;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,NeighbourEdges);
	for (vtkIdType i=0;i<NumberOfEdges;i++)
	{
		Edge &edge = this->Edges[ NeighbourEdges[i] ];
		vtkIdType v2=edge.Vertex1;
		if (v2==v1)	Output->InsertNextId(edge.Vertex2);
		else		Output->InsertNextId(v2);
	}
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::GetNeighbours(vtkIdList *Input,vtkIdList *Output)
{
	Output->Reset();
	vtkIdType NumberOfEdges,*NeighbourEdges;

	for (vtkIdType i=0;i<Input->GetNumberOfIds();i++)
	{
		vtkIdType v1 = Input->GetId(i);
		Output->InsertUniqueId(v1);
		this->GetVertexNeighbourEdges(v1,NumberOfEdges,NeighbourEdges);
		for (vtkIdType j=0;j<NumberOfEdges;j++)
		{
			Edge &e = this->Edges[ NeighbourEdges[j] ];
			vtkIdType v2 = e.Vertex1;
			if (v2==v1) Output->InsertUniqueId(e.Vertex2);
			else		Output->InsertUniqueId(v2);
		}
	}
}


vtkIdType vtkSurfaceBase::IsEdgeBetweenFaces(const vtkIdType &f1, const vtkIdType &f2)
{
	vtkIdType *Vertices1;
	vtkIdType NumberOfVertices1;
	this->GetFaceVertices(f1,NumberOfVertices1,Vertices1);
	vtkIdType *Vertices2;
	vtkIdType NumberOfVertices2;
	this->GetFaceVertices(f2,NumberOfVertices2,Vertices2);
	
	vtkIdType Vertices[2];
	
	int Count=0;
	for (vtkIdType i=0;i<NumberOfVertices1;i++)
	{
		vtkIdType v1=Vertices1[i];
		for (vtkIdType j=0;j<NumberOfVertices2;j++)
		{
			vtkIdType v2=Vertices2[j];
			if (v1==v2)
			{
				Vertices[Count++]=v1;
				break;
			}
		}
		if (Count==2)
			return (this->IsEdge(Vertices[0],Vertices[1]));		
	}
	return (-1);
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::InsertEdgeInRing(const vtkIdType &e1,const vtkIdType &v1)
{

	VertexRing &r=this->VerticesAttributes[v1];

	// First test if the edge is not already in ring
	int NumberOfEdges = r.size();
	for (vtkIdType i=0;i<NumberOfEdges;i++)
		if (r[i]==e1) return;

	r.push_back(e1);
}

void vtkSurfaceBase::DeleteEdgeInRing(const vtkIdType &e1, const vtkIdType &v1)
{
	VertexRing &r = this->VerticesAttributes[v1];
	int nEdges = r.size();
	vtkIdType i;
	for (i=0;i<nEdges;i++)
		if (r[i]==e1) {
			r[i]=r[nEdges-1];
			r.pop_back();
			return;
		}

}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::InsertFaceInRing(const vtkIdType &Face, const vtkIdType &e)
{
	vtkIdType F1,F2;
	this->GetEdgeFaces(e,F1,F2);
	Edge &edge = this->Edges[ e ];
	if (F1==-1)
	{
		edge.Poly1 = Face;
		return;
	}
	if (F1==Face) return;

	if (F2==-1)
	{
		edge.Poly2 = Face;
		return;
	}
	if (F2==Face) return;

	vtkIdList *FList=edge.NonManifoldFaces;
	if (!FList)
	{
		FList=vtkIdList::New();
		edge.NonManifoldFaces=FList;
	}
	FList->InsertUniqueId(Face);
	return;
}

void vtkSurfaceBase::DeleteFaceInRing(const vtkIdType &Face, const vtkIdType &e)
{
	Edge &edge = this->Edges[e ];
	vtkIdList *FList=edge.NonManifoldFaces;
	if (FList==0)
	{
		if (edge.Poly1==Face)
			edge.Poly1=edge.Poly2;
		edge.Poly2=-1;
		return;
	}
	vtkIdType Index=FList->IsId(Face);
	if (Index>=0)
		FList->DeleteId(Face);        
	else
	{
		int LastPosition=FList->GetNumberOfIds()-1;
		vtkIdType F3=FList->GetId(LastPosition);
		FList->DeleteId(F3);
		if (edge.Poly1==Face)
			edge.Poly1=F3;
		else
			edge.Poly2=F3;
	}
	if (FList->GetNumberOfIds()==0)
	{
		FList->Delete();
		edge.NonManifoldFaces=0;
	}
}

vtkIdType vtkSurfaceBase::AddVertex(double x, double y, double z)
{
	vtkIdType v1;
	if (this->VerticesGarbage.empty()) {
		v1=this->GetPoints()->InsertNextPoint(x,y,z);
		this->AllocateVerticesAttributes(v1 + 1);
	} else {
		v1=this->VerticesGarbage.front();
		this->VerticesGarbage.pop();
		this->VerticesAttributes[v1].resize(0);
		this->GetPoints()->SetPoint(v1,x,y,z);
	}
	this->GetPoints()->Modified();
	this->ActiveVertices->SetValue(v1,1);
	return v1;

}

// ****************************************************************
// ****************************************************************
vtkIdType vtkSurfaceBase::AddEdge(const vtkIdType &v1, const vtkIdType &v2,const vtkIdType &f1)
{
	if (v1==v2)
	{
		cout<<"Error : creation of a self-loop for vertex "<<v1<<endl;
		return (-1);
	}

	vtkIdType edge=this->IsEdge(v1,v2);

	if (edge>=0)
	{
		Edge &e = this->Edges[edge];
		if (e.Poly1<0)
		{
			e.Poly1=f1;
			return (edge);
		}

		if (e.Poly2>=0)
		{
			vtkIdList *LIST=e.NonManifoldFaces;
			if (!LIST)
			{
				LIST=vtkIdList::New();
				e.NonManifoldFaces=LIST;
			}
			LIST->InsertNextId(f1);
			return (edge);
		}
		else
		{
			e.Poly2=f1;
			return (edge);
		}
	}

	if (this->EdgesGarbage.empty()) {
		edge = this->Edges.size();
		this->Edges.emplace_back();
	} else {
		edge=this->EdgesGarbage.front();
		this->EdgesGarbage.pop();
	}

	Edge &e = this->Edges[edge];
	e.Vertex1=v1;
	e.Vertex2=v2;
	e.Poly1=f1;
	e.Poly2=-1;
	e.NonManifoldFaces=0;
	e.Active=true;
	this->InsertEdgeInRing(edge,v1);
	this->InsertEdgeInRing(edge,v2);
	return (edge);
}

vtkIdType vtkSurfaceBase::AddFace(const vtkIdType& v1,const vtkIdType& v2,const
vtkIdType& v3)
{
	int Number=3;
	vtkIdType Vertices[3];
	Vertices[0]=v1;
	Vertices[1]=v2;
	Vertices[2]=v3;
	return this->AddPolygon(Number,Vertices);
}

// ****************************************************************
// ****************************************************************
vtkIdType vtkSurfaceBase::AddPolygon(int NumberOfVertices,vtkIdType *Vertices)
{
	vtkIdType face;
	int Type=0;
	
	// Check if the vertices of the polygon have been created
	
	for (int i=0;i<NumberOfVertices;i++)
	{
		if (Vertices[i]>=this->GetNumberOfPoints())
			Type=99;		
	}
	if (Type)
	{
		cout<<"ERROR : attempt to create the Cell with vertices : ";
		for (int i=0;i<NumberOfVertices;i++)
			cout<<Vertices[i]<<" ";
		cout<<endl;
		cout<<"But Only "<<this->GetNumberOfPoints()<<" vertices have been created. Exiting program"<<endl;
		exit(1);
	}

	
	switch (NumberOfVertices)
	{
	case 3:
		Type=VTK_TRIANGLE;
		break;
	case 4:
		Type=VTK_QUAD;
		break;

	default:
		Type=VTK_POLYGON;
		break;
	}

	if (this->CellsGarbage[NumberOfVertices].empty())
	{
#if ( (VTK_MAJOR_VERSION < 9))
		face=this->Polys->InsertNextCell(NumberOfVertices);
		for (int i=0;i<NumberOfVertices;i++)
			this->Polys->InsertCellPoint(Vertices[i]);
		int Loc=this->Polys->GetInsertLocation(NumberOfVertices);
		this->Cells->InsertNextCell(Type,Loc);
#else
		face=this->Polys->InsertNextCell(NumberOfVertices, Vertices);
#endif
		this->AllocatePolygonsAttributes(face + 1);
	}
	else
	{
		face=this->CellsGarbage[NumberOfVertices].front();
		this->CellsGarbage[NumberOfVertices].pop();
#if ( (VTK_MAJOR_VERSION < 9))
		int loc = this->Cells->GetCellLocation(face);
		this->Cells->InsertCell(face,Type,loc);
		this->Polys->ReplaceCell(loc,NumberOfVertices,Vertices);
#else
		this->Polys->ReplaceCellAtId(face,NumberOfVertices,Vertices);
#endif
	}

	for (int i=0;i<NumberOfVertices;i++)
	{
		this->AddEdge(Vertices[i],Vertices[(i+1)%NumberOfVertices],face);
	}

	if ((this->FirstTime)&&(face==0))
	{
		this->VisitedPolygons->SetValue(face,1);
		this->FirstTime=false;
	}
	else
		this->VisitedPolygons->SetValue(face,0);

	this->ActivePolygons->SetValue(face,1);
	if (this->OrientedSurface)
	{
		this->ConquerOrientationFromFace(face);
	}

	this->Polys->Modified();
	this->Modified();
	return (face);
}

// ****************************************************************
// ****************************************************************


void vtkSurfaceBase::AllocateVerticesAttributes(int NumberOfVertices)
{
	if ( NumberOfVertices <= this->VerticesAttributes.size() ) return;
	this->VerticesAttributes.resize(NumberOfVertices) ;
	this->ActiveVertices->Resize(NumberOfVertices);
	this->VerticesAttributes.reserve( NumberOfVertices );
}

void vtkSurfaceBase::AllocatePolygonsAttributes(int NumberOfPolygons)
{
	if (this->VisitedPolygons->GetSize()>=NumberOfPolygons) return;
	this->VisitedPolygons->Resize(NumberOfPolygons);
	this->ActivePolygons->Resize(NumberOfPolygons);

}

// ****************************************************************
// ****************************************************************
// cette fonction initialise l'objet a partir d'un autre object de meme type
// Alloue la place et l'objet pour stocker le meme nombre de points de
// triangles et d'arretes , MAIS NE COPIE PAS LES DONNEES
void vtkSurfaceBase::Init(vtkSurfaceBase *mesh)
{

	// recupere le nombre de triangles, points et arretes
	int numPoints=mesh->GetNumberOfPoints();
	int numFaces=mesh->GetNumberOfCells();
	int numEdges=mesh->GetNumberOfEdges();

	// delegue le travail
	this->Init(numPoints,numFaces,numEdges);

}


// ****************************************************************
// ****************************************************************
// cette fonction initialise l'objet a partir d'un nombre de points,
// de triangles et d'arretes , renvoie un objet vide
void vtkSurfaceBase::Init(int numPoints, int numFaces, int numEdges)
{

	vtkCellArray *CellsArray1;
	CellsArray1 = vtkCellArray::New();
	#if ( (VTK_MAJOR_VERSION < 9))
	CellsArray1->Allocate(4*numFaces,numFaces);
	#else
	CellsArray1->Allocate(3*numFaces,numFaces);
	#endif
	this->SetPolys(CellsArray1);
	CellsArray1->Delete();
	#if ( (VTK_MAJOR_VERSION >= 9))
	CellsArray1->Use64BitStorage();
	#else
	if (this->Cells)
		this->Cells->Delete();

	this->Cells = vtkCellTypes::New();
	this->Cells->Allocate(numFaces,numFaces);
	#endif

	// create and allocate memory for points
	vtkPoints *Points1=vtkPoints::New();
	Points1->Allocate(numPoints);
	this->SetPoints(Points1);
	Points1->Delete();

	// create and allocate memory for all the vtkSurfaceBase specific tables
	this->AllocateVerticesAttributes(numPoints);
	this->AllocatePolygonsAttributes(numFaces);
	this->Edges.reserve(numEdges);

}

// ****************************************************************
// ****************************************************************
// fonction CreateFromPolyData
// cette fonction initialise l'objet a partir d'un polydata
// copie les donnees (points et triangles)
// et cree le tableau d'arretes
void vtkSurfaceBase::CreateFromPolyData(vtkPolyData *input)
{
	vtkIdType i,j,v1,v2;
	vtkIdType NumberOfVertices;
	vtkIdType *Vertices;

	// just copy the polydata in input
	this->ShallowCopy(input);

	// Delete the cells that are not polygons
	vtkCellArray *VerticesCells=this->GetVerts();
	vtkCellArray *LinesCells=this->GetLines();
	vtkCellArray *StripsCells=this->GetStrips();

	if ((VerticesCells!=0)||(LinesCells!=0)||(StripsCells!=0))
	{
		this->SetVerts(0);
		this->SetLines(0);
		this->SetStrips(0);
		this->Modified();
		this->BuildCells();
	}

	vtkIdType numPoints=this->GetNumberOfPoints();
	vtkIdType numFaces=this->GetNumberOfCells();

	this->AllocateVerticesAttributes(numPoints);
	this->AllocatePolygonsAttributes(numFaces);

	for (i=0;i<this->GetNumberOfCells();i++)
	{
		bool ActiveFace=false;
		this->GetFaceVertices(i,NumberOfVertices,Vertices);
		if (NumberOfVertices>1)
		{
			// test whether the face is an active one (at least its first two vertices should be different)
			if (Vertices[0]!=Vertices[1])
			{
				ActiveFace=true;
				for (j=0;j<NumberOfVertices;j++)
				{
					v1=Vertices[j];
					v2=Vertices[(j+1)%NumberOfVertices];
					this->AddEdge(v1,v2,i);
				}
			}
		}
		if (ActiveFace)
			this->ActivePolygons->SetValue(i,1);
		else
		{
			// The face is not used. Let's push it in the garbage collector
			this->ActivePolygons->SetValue(i,0);
			this->CellsGarbage[NumberOfVertices].push(i);			
		}
	}

	if (this->OrientedSurface)
	{
		this->CheckNormals();
	}
}

vtkSurfaceBase::vtkSurfaceBase()
{
	this->FirstTime=true;
	this->SetOrientationOn();

	// vertices attributes
	this->ActiveVertices = vtkBitArray::New();

	// faces attributes 
	this->VisitedPolygons = vtkBitArray::New();
	this->ActivePolygons = vtkBitArray::New();

	this->CleanEdges=1;
	this->CleanVertices=0;
	
	this->Init(50,100,150);
}

// ****************************************************************
// ****************************************************************
vtkSurfaceBase::~vtkSurfaceBase() {

	for ( auto it = this->Edges.begin(); it != this->Edges.end(); it++)
		if ( it->NonManifoldFaces ) it->NonManifoldFaces->Delete();

	this->VisitedPolygons->Delete();
	this->ActivePolygons->Delete();
	this->ActiveVertices->Delete();
}
