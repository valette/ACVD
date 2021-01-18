/*=========================================================================

  Program:   Mailleur 3D multi-resolution
  Module:    vtkSurface.h
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

#ifndef __vtkSurface_h
#define __vtkSurface_h
#include <sstream>
#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkCell.h>
#include <vtkEdgeTable.h>
#include <vtkCommand.h>
#include <vtkIdListCollection.h>
#include <vtkPolyDataWriter.h>
#include <vtkPLYWriter.h>
#include "vtkSurfaceBase.h"


/**
 * This class provides easy access and update for polygonal meshes
 * It is actually derived from the vtkSurfaceBase which handles all the core functionnalities
 * vtkSurface only adds general purpose methods: load a mesh from a file, compute its properties,
 * subdivide etc....
 * See vtkSurfaceBase.h for core functionnalities
 *
 */
class VTK_EXPORT vtkSurface : public vtkSurfaceBase
{
public:

	/// returns a vtkSurface with no empty memory slots
	vtkSurface *CleanMemory();

	/// Adds uniform noise to the model. The noise distribution is uniform, and its strength
	/// depends on Magnitude*BoundingBoxDiagonalLength
	void AddNoise(double Magnitude);

	/// Returns the points bounding box diagonal length
	double GetBoundingBoxDiagonalLength();

	/// Recusrsively splits the edges longer than Ratio*AverageLength
	void SplitLongEdges (double Ratio);

	/// returns a mesh which is a linear subdivision of this
	vtkSurface *Subdivide(vtkIntArray *Parent1=0, vtkIntArray *Parent2=0);

	/// subdivide the mesh (in place)
	void SubdivideInPlace(vtkIntArray *Parent1=0, vtkIntArray *Parent2=0);

	/// copy to stream the connectivity and geometric properties of the mesh
	void GetMeshProperties(std::stringstream &stream);

	/// same as GetMeshProperties, but the properties are displayed on screen
	void DisplayMeshProperties();
	
	/// Switch the whole mesh orientation if its normals are not directed outwards.
	/// This is done by computing the signed volume of the mesh
	void EnsureOutwardsNormals();

	/// Computes the normal of the given triangle
	void GetTriangleNormal(vtkIdType Triangle, double *Normal);

	/// Computes the normal at a vertex (works only for triangular meshes)
	void GetVertexNormal(vtkIdType Vertex, double *Normal);
	
	/// returns the legth of the input edge;
	double GetEdgeLength(vtkIdType Edge);

	/// returns the normals of the mesh triangles
	vtkDoubleArray *GetTrianglesNormals();
	void DeleteTrianglesNormals();

	// returns the area of a given cell and its barycenter
	void GetCellMassProperties(vtkIdType CellId, double &Area, double *Baricenter);

	/// Computes the Triangles Areas
	vtkDoubleArray * GetTrianglesAreas();
	void DeleteTrianglesAreas();

	/// Returns the area of the Face
	double GetFaceArea(vtkIdType Face);

	/// returns a polydata of the mesh edges only
	vtkPolyData *GetEdgesPolyData();

	/// returns a polydata of the mesh vertices only
	vtkPolyData *GetVerticesPolyData();

	/// rescale the points so that the maximum coordinates range is equal to 1
	void RescaleCoordinates();

	/// Returns the euclidian distance between the vertices V1 and V2
	double GetDistanceBetweenVertices(vtkIdType V1,vtkIdType V2);

	/// Builds a 1 connected component mesh from the mesh in two steps:
	/// First: merge the points which are closer than the given relative tolerance
	/// Second : add fictive edges between unconnected components
	vtkSurface* CleanConnectivity(double tolerance);

	/// Computes the area of Vertex (it will actually compute the area of the
	/// Surrounding cells divided by their respective number of vertices)
	double GetVertexArea(vtkIdType Vertex);

	/// Computes the Vertices Areas with respect to their surrounding cells
	vtkDoubleArray* GetVerticesAreas();
	void DeleteVerticesAreas();

	/// Computes the Mesh Edges Lengths in an array
	vtkDoubleArray* GetEdgeLengths();
	void DeleteEdgeLengths();

	/// Computes the connected components of the mesh
	vtkIdListCollection* GetConnectedComponents();
	void DeleteConnectedComponents();

	/// returns a vtkSurface made of the biggest connected component
	vtkSurface *GetBiggestConnectedComponent();

	/// returns a vtkSurface made of the n biggest connected components
	vtkSurface *GetBiggestConnectedComponents( int numberOfComponents );

   /// Computes the mesh minimal angle, average minimal angle, minimal triangle quality and average triangle quality
   /// usefull to give an objective quality criterion of the mesh
   void ComputeTrianglesStatistics(double &Amin, double & Aav,double &Qmin, double &Qav, double &P30);

   /// Computes the histogram of triangles quality and writes it in a file
   void ComputeQualityHistogram(const char *FileName);

   /// Creates the vtkSurface object by reading a mesh file.
   /// Supported file types : .wrl, .vtk and .ply
   void CreateFromFile (const char *FileName);

   /// Saves the mesh
   /// Supported file types : .vtk and .ply
   void WriteToFile (const char *FileName);

   /// The Constructor vtkSurface::New();
   static vtkSurface *New();

   vtkTypeMacro(vtkSurface,vtkSurfaceBase);

   /// Create a similar type object.
   vtkDataObject *MakeObject()
     {return vtkSurface::New();};

   /// Quantizes the coordinates of the mesh to integers
   /// (dynamic range : q bits)
   void QuantizeCoordinates(int q);

   /// Quantizes the coordinates of the mesh to integers
   /// (quantization parameters will be the same as Input
   void QuantizeCoordinatesLike(vtkSurface *Mesh);

     /// Quantizes the coordinates of the mesh to integers
   /// (quantization parameters must be given
   void QuantizeCoordinates(double Factor, double Tx, double Ty, double Tz);


   /// Writes the mesh to a .iv file
   void WriteInventor(const char *filename);

   /// Writes the mesh to a .smf file
   void WriteSMF(const char *filename);

   /// Computes the list of sharp vertices
   void ComputeSharpVertices(double treshold);

   /// Returns 1 if v is a sharp vertex
   /// Otherwise Returns 0.
   int IsSharpVertex(vtkIdType v)
     {return (this->SharpVertices->GetValue(v));};

   /// Cancels the computation of the sharp vertices list
   void DeleteSharpVertices() ;


   /// returns the scaling factors used for quantization
   void GetScalingFactors (double &Factor, double &Tx, double &Ty, double &Tz)
     {
       Factor=this->Factor;
       Tx=this->Tx;
       Ty=this->Ty;
       Tz=this->Tz;
     };

   /// sets the scaling factors used for quantization
   void SetScalingFactors (double Factor, double Tx, double Ty, double Tz)
     {
       this->Factor=Factor;
       this->Tx=Tx;
       this->Ty=Ty;
       this->Tz=Tz;
     };

   /// sets the scaling factors used for quantization
   void SetScalingFactors (vtkSurface *Mesh)
     {
	double Factor,Tx,Ty,Tz;
	Mesh->GetScalingFactors(Factor,Tx,Ty,Tz);
       this->Factor=Factor;
       this->Tx=Tx;
       this->Ty=Ty;
       this->Tz=Tz;
     };

   void UnQuantizeScalingFactors ();
   void UnQuantizeCoordinates();

	void PrintVerticesCoordinates();
	void PrintConnectivity();
	void SaveConnectivity(const char * FileName);

    void WriteMeshStatisticsFile( const char* AreaFileName, const char* QFileName,
                             const char* AngleMinFileName, const char* DegreeFileName );

protected:

	vtkSurface();
	~vtkSurface();

private:

   /// Attribute on the vertices:
   vtkIntArray *SharpVertices;

   vtkDoubleArray *VerticesAreas;
   vtkDoubleArray *TrianglesAreas;
   vtkDoubleArray *TrianglesNormals;

	vtkDoubleArray *EdgeLengths;


	/// parameters used for quantization
   double Tx;
   double Ty;
   double Tz;
   double Factor;

	/// A Collection of IdLists containing Ids of each connected component
	vtkIdListCollection* ConnectedComponents;


};


#endif
