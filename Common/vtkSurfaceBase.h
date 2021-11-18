/*=========================================================================

  Program:   vtkSurfaceBase
  Module:    vtkSurface
  Language:  C++
  Date:      2002/05
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

#ifndef __vtkSurfaceBase_h
#define __vtkSurfaceBase_h

#include <map>
#include <queue>

#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkBitArray.h>
#include <vtkPolyData.h> 
#include <vtkCell.h>
#include <vtkCommand.h>
#include <vtkCellArray.h>

/**
 *  An efficient class for 3D triangular mesh processing.
 *  The vtkPolyData data structure wasn't designed for manipulating
 *  triangular meshes (but only to visualize them). It misses a true notion
 *  of edge with efficient accesses to adjacent vertices, adjacent faces,
 *  and more generaly the traversal of the graph of edges (like in half-edge
 *  data structure).
 *  The vtkSurfaceBase Class is an attempt to represent a low-level
 *  (because the internal representation is based on arrays and not
 *  pointers) edge-oriented 3D triangular mesh within the
 *  vtk hierarchy, and enables constant-time connectivity queries,
 *  edges labeling and easy update functionnalities.
 */

class VTK_EXPORT vtkSurfaceBase : public vtkPolyData
{

public:

	/// The Constructor vtkSurfaceBase::New();
	static vtkSurfaceBase *New();
	vtkTypeMacro(vtkSurfaceBase,vtkPolyData);
	
	/// Create a similar type object.
	vtkDataObject *MakeObject()
	{return vtkSurfaceBase::New();};
	
	/// Embeds a vtkPolyData into a vtkSurfaceBase object.
	void CreateFromPolyData(vtkPolyData *input);
	
	/// Set memory allocation in the vtkSurfaceBase object for further Cells
	/// and points insertion (this is not mandatory, and maybe useless)
	void Init (int numPoints, int numFaces, int numEdges);
	
	/// Set memory allocation in the vtkSurfaceBase object, which is assumed to
	/// be as large as the *mesh vtkSurfaceBase
	void Init (vtkSurfaceBase *mesh);

	/// Recover extra memory
	virtual void SQueeze();

	/// adds a vertex with coordinates (x,y,z) to the vtkSurface and returns its Id
	vtkIdType AddVertex(double x, double y, double z);

	/// adds a vertex with coordinates (x,y,z) to the vtkSurface and returns its Id
	vtkIdType AddVertex(double *P);

	/// adds an edge (v1,v2) to the vtkSurface and returns its Id
	vtkIdType AddEdge(const vtkIdType& v1,const vtkIdType& v2);

	/// Adds a face (ve1,ve2,ve3) to the vtkSurfaceBase Object and returns its Id
	vtkIdType AddFace(const vtkIdType& ve1, const vtkIdType& ve2, const vtkIdType& ve3);

	/// Adds a face (ve1,ve2,ve3) to the vtkSurfaceBase Object and returns its Id
	vtkIdType AddPolygon(int NumberOfVertices, vtkIdType*Vertices);

	/// Flips the edge. If the flip was performed, return -1. Otherwise, the edge already existing is returned.
	vtkIdType FlipEdge(vtkIdType edge);
	
	/// Flips the edge. If the resulting edge already exists, it is flipped (an so on). Returns -1 in case of success
	vtkIdType FlipEdgeSure(vtkIdType edge);

	/// Bisects Edge by creating a new vertex at its midpoints and two new faces.
	vtkIdType BisectEdge(vtkIdType e);

	/// Call this function if you want the surface to be oriented.
	/// By default, the surface is oriented.
	void SetOrientationOn();
	
	/// Call this function if you do not want the surface to be oriented.
	/// By default, the surface is oriented.
	void SetOrientationOff();

	/// Reorders Cells description so that each cell is well oriented.
	/// Usefull for wireframe rendering.
	/// Note : this function us usefull ONLY if you haven't called
	///        SetOrientationOn() before.
	void CheckNormals();

	/// Deletes the face f1, possibly its adjacent edges and vertices if they are free
	/// and is CleanVertices and CleanEdges are set to 1
	void DeleteFace(vtkIdType f1);

	/// Deletes the free edge e1 and possibly its vertices if they are free.
	/// free means not connected to any higher order element
	void DeleteEdge(vtkIdType e1);

	/// Deletes the FREE vertex v1
	/// free means not connected to any higher order element
	void DeleteVertex(vtkIdType v1);

	/// Sets On/Off the ability to clean vertices.  Default Value : 1
	void SetCleanVertices(int C) {this->CleanVertices=C;};

	/// returns the value of CleanVertices. Default Value : 1
	int GetCleanVerticesState(){return this->CleanVertices;};

	/// returns the value of CleanEdges. Default Value : 1
	int GetCleanEdgesState(){return this->CleanEdges;};

	/// Sets On/Off the ability to clean Edges. Default Value : 1
	void SetCleanEdges(int C) {this->CleanEdges=C;};

	/// Merge the two input vertices v1 and v2. This is a general general merge operator
	/// which can be used as edge collapse operator. (NOT FULLY TESTED)
	void MergeVertices(vtkIdType v1, vtkIdType v2);

	/// Changes one vertex of a given face
	void ChangeFaceVertex(vtkIdType Face, vtkIdType OldVertex, vtkIdType NewVertex);

	/// Returns the number of edges in the vtkSurfaceBase Object
	int GetNumberOfEdges() {return this->Edges.size();};
	
	/// Returns v1 and v2 as the vertices bounding the edge 
	void GetEdgeVertices(const vtkIdType& edge, vtkIdType &v1, vtkIdType &v2);

	// returns the number of faces adjacvent to the input edge
	int GetEdgeNumberOfAdjacentFaces(const vtkIdType &e);
	
	/// Returns the vertices v1, v2 and v3 bounding the face
	void GetFaceVertices(const vtkIdType& face,
		vtkIdType &v1, vtkIdType &v2, vtkIdType &v3);

	/// Returns the number of vertices and a pointer to them for the face 
	void GetFaceVertices(const vtkIdType& face, vtkIdType &NumberOfVertices,
vtkIdType* &Vertices);
	
	/// Returns the coordinates of Point
	void GetPointCoordinates(vtkIdType Point, double *x);
	
	/// Returns the valence of the Point V1 i.e. its number of adjacent Points
	int GetValence(const vtkIdType& v1);
	
	/// Returns the number of boundaries at Point v1
	int GetNumberOfBoundaries(const vtkIdType &v1);
	
	/// Conquer from a triangle f1 through the edge bounded by v1 and v2. 
	/// f2 is the conquered triangle and v3 its third vertex 
	/// (v1 and v2 are the first and second vertices)
	void Conquer(const vtkIdType& f1,const vtkIdType& v1, const vtkIdType& v2,
		vtkIdType &f2, vtkIdType &v3);
	
	/// Returns f1 and f2 as the faces adjacent to the edge.
	/// WARNING : if the edge is a boundary edge, f2=-1 (it has only
	///           1 adjacent face)
	void GetEdgeFaces(const vtkIdType& e1, vtkIdType &f1, vtkIdType &f2);

	/// Returns the faces adjacent to the edge as an IdList
	void GetEdgeFaces(vtkIdType e1, vtkIdList *Flist);
	
	/// Returns true if the edge is a Manifold one (2 adjacent polygons)
	bool IsEdgeManifold(const vtkIdType& e);

	/// returns true if the vertex is manifold.
    bool IsVertexManifold( const vtkIdType& iV);
	
	/// Returns the List of Points adjacent to the point v1	
	void GetVertexNeighbours(vtkIdType v1, vtkIdList *Output);
	
	/// Returns the List of Faces adjacent to the Face	
	void GetFaceNeighbours(vtkIdType Face,vtkIdList *FList);

	/// Returns the List of Faces adjacent to the Face	
	void GetFaceNeighbours(vtkIdType Face,vtkIdListCollection *FList);
	
	/// Returns the list of edges adjacent to the point v1
	void GetVertexNeighbourEdges(vtkIdType v1, vtkIdList *Output);

	/// Returns the pointer to the edges in the ring of v1. Faster than the method with vtkIdList
	/// Be carefull not to modify the pointed values, it would break the structure!
	/// Note that when the mesh is modified, the data pointer might get invalidated
	void GetVertexNeighbourEdges(const vtkIdType& v1, vtkIdType &NumberOfEdges,
vtkIdType* &Edges);

	/// Returns the list of faces adjacent to the point v1
	void GetVertexNeighbourFaces(const vtkIdType &v1, vtkIdList *Output);
	
	/// Returns the list of vertices adjacent to the list of vertices *Input.
	/// Usefull for constructing n-rings
	void GetNeighbours(vtkIdList *Input, vtkIdList *Output);
	
	/// Returns the Id of the face having v1, v2 and v3 for vertices.
	/// Returns -1 if the face doesn't exist.
	vtkIdType IsFace(const vtkIdType &v1, const vtkIdType &v2, const vtkIdType &v3);
	
	/// Returns the Id of the edge (v1,v2).
	/// Returns -1 if the edge doesn't exist.
	vtkIdType IsEdge(const vtkIdType &v1, const vtkIdType &v2);
	
	/// Returns the Id of the edge between the faces f1 and f2
	/// Returns -1 if the edge doesn't exist
	vtkIdType IsEdgeBetweenFaces(const vtkIdType &f1, const vtkIdType &f2);
	
	/// Returns the third Point of the face f1 (the two first Points are v1 and v2) 
	vtkIdType GetThirdPoint(const vtkIdType& f1,const vtkIdType& v1,const vtkIdType& v2);
	
	/// Sets the coordinates of Point
	void SetPointCoordinates(const vtkIdType &Point, double *x);

	/// Returns the first edge in the ring of v1
	vtkIdType GetFirstEdge(const vtkIdType& v1);
	
	/// Returns a boundary edge in the ring of v1. If there is no boundary edge, any given edge is returned
	vtkIdType GetBoundaryEdge(const vtkIdType& v1);

	/// Returns the entropy of the valences of the mesh
	double GetValenceEntropy();
	
	/// writes a text file containing an histogram of the vertices valence distribution
	void GetValenceTab(char *FileName);

	/// switches the cells orientation (usefull when the mesh is displayed all black...)
	void SwitchOrientation();
	
	/// returns 1 if Vertex is actually used to store a polygon (not deleted). Returns 0 otherwise
	int IsVertexActive(const vtkIdType &Vertex) {return (this->ActiveVertices->GetValue(Vertex));};

	/// returns 1 if Face is actually used to store a polygon (not deleted). Returns 0 otherwise
	int IsFaceActive(const vtkIdType &Face) {return (this->ActivePolygons->GetValue(Face));};

	/// returns 1 if Edge is actually used to store an edge (not deleted). Returns 0 otherwise	
	int IsEdgeActive(const vtkIdType &e) {return this->Edges[e].Active;};
	
	/// Checks the integrity of the structure. Returns true if the structure is OK
	bool CheckStructure();

	/// Outputs the underlying data in readable format (not suited or large meshes)
	void DisplayInternals( bool exitOnError = false );

protected:

	/// the constructor
	vtkSurfaceBase(); 

	/// the desctructor
	~vtkSurfaceBase();

private:

	class Edge {
		public:
		/// Edges attributes:
		/// - one of the face adjacent to this edge
		vtkIdType Poly1;
		/// - the other face adjacent to this edge
		vtkIdType Poly2;
		/// - one of the vertex adjacent to this edge
		vtkIdType Vertex1;
		/// - the other vertex adjacent to this edge
		vtkIdType Vertex2;

		///  lists of other faces adjacent to the edges.
		///  if no more than two faces are adjacent to the edge, then EdgesNonManifoldFaces=0
		vtkIdList *NonManifoldFaces;
		// This array determines whether an edge slot is used or not
		bool Active;
		Edge() {
			Poly1 = Poly2 = Vertex1 = Vertex2 = -1;
			NonManifoldFaces = 0;
			Active = false;
		}
	};

	/// adds an edge (v1,v2) to the vtkSurface, with possibly an adjacent face f1
	vtkIdType AddEdge(const vtkIdType &v1, const vtkIdType &v2, const vtkIdType &f1);

	/// Replaces the definition of f1 by the face formed by v1,v2 and v3.
	/// WARNING : f1 must already exist.
	void SetFace(const vtkIdType& f1,const vtkIdType& v1,const vtkIdType& v2,const vtkIdType& v3);

	/// inserts the edge e1 in the ring of v1
	void InsertEdgeInRing(const vtkIdType &e1, const vtkIdType &v1);
	/// removes the edge e1 from the ring of v1
	void DeleteEdgeInRing(const vtkIdType &e1, const vtkIdType &v1);

	///  inserts Face in the ring of Edge
	void InsertFaceInRing(const vtkIdType &Face, const vtkIdType &e);
	///  removes Face from the ring of Edge
	void DeleteFaceInRing(const vtkIdType &Face, const vtkIdType &e);

	/// deletes the edge if it is adjacent to no polygon
	/// uneffective if CleanEdges is set to 0
	void CleanEdge(const vtkIdType &e);

	/// deletes the vertex if it is not connected
	/// uneffective if CleanVertices is set to 0
	void CleanVertex(const vtkIdType &Vertex);

	///=true if we keep orientation of the surface
	bool OrientedSurface;
	void AllocateVerticesAttributes(int NumberOfPoints);
	void AllocatePolygonsAttributes(int NumberOfPolygons);

	/// VertexRing contains all adjacent edges to a given vertex
	typedef std::vector< vtkIdType > VertexRing;
	std::vector <VertexRing> VerticesAttributes;
	/// This array determines whether a vertex slot is used or not
	vtkBitArray *ActiveVertices;

	std::vector <Edge> Edges;

	/// Array determining the polygons which were already visited by the consistent-orientation keeping method
	vtkBitArray *VisitedPolygons;
	/// This array determines whether a polygon slot is used or not
	vtkBitArray *ActivePolygons;

	/// this flag is used only once, to set the face 0 as visited (and not after, when it is removed)
	bool FirstTime;

	/// the method which keeps the orientation of the faces consistent. This is called whenever a new face is created
	void ConquerOrientationFromFace(vtkIdType Face);

	/// Parameter to control the Vertices cleaning step when deleting faces
	int CleanVertices;

	/// Parameter to control the Edges cleaning step when deleting faces
	int CleanEdges;

	/// Garbage collector for deleted vertices
	std::queue<vtkIdType> VerticesGarbage;
	/// Garbage collector for deleted cells
	std::map<int, std::queue<int>> CellsGarbage;
	/// Garbage collector for deleted edges
	std::queue<int> EdgesGarbage;
};

inline void vtkSurfaceBase::SetFace(const vtkIdType& f1,
									const vtkIdType& v1,const vtkIdType& v2,const vtkIdType& v3)
{
	vtkIdType n, *vertices;
	this->GetFaceVertices( f1, n, vertices );
	vertices[0]=v1;
	vertices[1]=v2;
	vertices[2]=v3;
	this->Polys->Modified();
};
// ** METHODE GetThirdPoint
inline vtkIdType vtkSurfaceBase::GetThirdPoint(const vtkIdType& f1,const vtkIdType& v1,const vtkIdType& v2)
{
	vtkIdType v3,v4,v5;
	this->GetFaceVertices(f1,v3,v4,v5);

	if ((v1!=v3)&&(v2!=v3))
	{
		return (v3);
	}
	else
	{
		if ((v1!=v4)&&(v2!=v4))
			return (v4);
		else
			return (v5);
	}
}

inline void vtkSurfaceBase::Conquer(const vtkIdType& f1,const vtkIdType& v1,const vtkIdType& v2,
									vtkIdType &f2, vtkIdType &v3)
{
	vtkIdType edge=this->IsEdge(v1,v2);
	if (edge<0)
	{
		f2=-1;
		v3=-1;
		return;
	}

	f2=this->Edges[edge].Poly2;
	if (f2==-1)
	{
		v3=-1;
		return;
	}

	if (f2==f1)
	{
		f2=this->Edges[edge].Poly1;
	}
	v3=this->GetThirdPoint(f2,v1,v2);
}

inline void vtkSurfaceBase::GetEdgeVertices(const vtkIdType& edge, vtkIdType &v1,
											vtkIdType &v2)
{
	Edge &e = this->Edges[ edge ];
	v1 = e.Vertex1;
	v2 = e.Vertex2;
}

inline void vtkSurfaceBase::GetFaceVertices(const vtkIdType &face, vtkIdType &v1, vtkIdType &v2, vtkIdType &v3)
{
	vtkIdType n;
	vtkIdType *pts;
	this->GetFaceVertices( face, n, pts );
	v1=pts[0];
	v2=pts[1];
	v3=pts[2];
}

#if ( (VTK_MAJOR_VERSION < 9))
inline void vtkSurfaceBase::GetFaceVertices(const vtkIdType& face,
vtkIdType &NumberOfVertices, vtkIdType* &Vertices)
{
	this->Polys->GetCell( this->Cells->GetCellLocation( face ), NumberOfVertices, Vertices );
}
#else

inline void vtkSurfaceBase::GetFaceVertices(const vtkIdType& face,
vtkIdType &NumberOfVertices, vtkIdType* &Vertices)
{
	auto connectivity = this->Polys->GetConnectivityArray64();
	auto offsets = this->Polys->GetOffsetsArray64();
	auto offset1 = offsets->GetValue( face );
	auto offset2 = offsets->GetValue( face + 1 );
	NumberOfVertices = offset2 - offset1;
	Vertices = static_cast < vtkIdType *> ( connectivity->GetPointer( offset1 ) );
}

# endif
inline void vtkSurfaceBase::GetVertexNeighbourEdges(const vtkIdType& v1,
vtkIdType &NumberOfEdges, vtkIdType* &Edges)
{
	VertexRing &r = this->VerticesAttributes[v1];
	NumberOfEdges = r.size();
	Edges=r.data();
}

/// Compute the valence (number of adjacent edges) of the given.
inline int vtkSurfaceBase::GetValence(const vtkIdType& v1)
{
	return this->VerticesAttributes[v1].size();
}

inline vtkIdType vtkSurfaceBase::GetBoundaryEdge(const vtkIdType& v1)
{
	vtkIdType *Edges,NumberOfEdges;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
	for (vtkIdType i=0;i<NumberOfEdges;i++)
	{
		vtkIdType e=Edges[i];
		Edge &edge=this->Edges[ e ];
		if ((edge.Poly1>=0)&&(edge.Poly2==-1))
			return (e);
	}
	return (*Edges);
}

inline vtkIdType vtkSurfaceBase::AddVertex(double *P)
{
	return this->AddVertex(P[0],P[1],P[2]);
}

inline vtkIdType vtkSurfaceBase::AddEdge(const vtkIdType& v1,const vtkIdType& v2)
{
	return this->AddEdge(v1,v2,-1);
}

inline void vtkSurfaceBase::GetPointCoordinates(vtkIdType Point, double *x)
{
	this->Points->GetPoint(Point,x);
}

inline void vtkSurfaceBase::GetEdgeFaces(const vtkIdType& e1, vtkIdType &f1, vtkIdType &f2) 
{
	Edge &e = Edges[ e1 ];
	f1=e.Poly1;
	f2=e.Poly2;
}

inline bool vtkSurfaceBase::IsEdgeManifold(const vtkIdType& e)
{
	Edge &edge = this->Edges[ e ];
	if ( ( edge.Poly2<0)||(edge.NonManifoldFaces!=0))
		return (false);
	else
		return (true);
}

inline void vtkSurfaceBase::SetPointCoordinates(const vtkIdType &Point, double *x)
{
	this->Points->SetPoint(Point,x);
	this->Points->Modified();
}

inline vtkIdType vtkSurfaceBase::GetFirstEdge(const vtkIdType&  v1)
{
	VertexRing &r = this->VerticesAttributes[v1];
	if(!r.size()) return(-1);
	return r[0];
}

inline vtkIdType vtkSurfaceBase::IsEdge(const vtkIdType &v1, const vtkIdType &v2)
{
	vtkIdType NumberOfEdges,*EdgesArray;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,EdgesArray);

	if (NumberOfEdges<=0) return (-1);
	for (NumberOfEdges--;NumberOfEdges!=-1;NumberOfEdges--)
	{
		vtkIdType e = EdgesArray[NumberOfEdges];
		Edge &edge=this->Edges[ e ];
		if ((v2==edge.Vertex1)||(v2==edge.Vertex2))
			return (e);
	}
	return (-1);
}

#endif

