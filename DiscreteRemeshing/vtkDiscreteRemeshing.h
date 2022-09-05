/***************************************************************************
                          vtkDiscreteRemeshing.h  -  description
                             -------------------
    begin                : jeu sep 25 2003
    copyright            : (C) 2003 by Sebastien Valette
    email                : 
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
 
#ifndef _VTKDISCRETEREMESHING_H_
#define _VTKDISCRETEREMESHING_H_

#include <vtkObjectFactory.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPLYWriter.h>
#include <vtkMetaImageReader.h>
#include <vtkMINCImageReader.h>
#include <vtkImageData.h>

#include "vtkSurface.h"
#include "RenderWindow.h"
#include "vtkSurfaceClustering.h"
#include "vtkCurvatureMeasure.h"
#include "vtkTag.h"

// A Class to process coarsening of vtkSurface PolyData. It derives from vtkSurfaceClustering
// This class adds meshing features to the clustering class, as proposed in : 
// " Approximated Centroidal Voronoi Diagrams for Uniform Polygonal Mesh Coarsening", Valette & Chassery, Eurographics 2004

template < class Metric > class vtkDiscreteRemeshing:public vtkSurfaceClustering <Metric >

{
    public:

	/// returns the coarsened model.
	vtkGetMacro( Output, vtkSurface* );

	// process the remeshing
	virtual void Remesh();

	/// for debugging: sets on/off the capability to skip the curvature computation and/or the clustering
	vtkSetMacro( FileLoadSaveOption, int )

	/// defines the Subsampling threshold. If the subsampling ratio is below this threshold,
	/// the mesh will be subdivided accordingly. Default value: 10
	vtkSetMacro( SubsamplingThreshold, int )

	/// Sets On/Off the Edges optimization scheme (still experimental)
	vtkSetMacro( EdgeOptimization, int )
	
	/// Sets On/Off the fix for meshes with boundaries. Default value: 1 (On)
	vtkSetMacro( BoundaryFixing, int )

	vtkSetMacro( ForceManifold, bool )
	vtkSetMacro( InputDensityFile, char* )
	vtkSetMacro( MaxCustomDensity, double )
	vtkSetMacro( MinCustomDensity, double )
	vtkSetMacro( CustomDensityMultiplicationFactor, double )

protected:

	vtkDiscreteRemeshing();
	~vtkDiscreteRemeshing();

	/// In this method, we compute the curvature indicator and adapt it to the possibly subdivided Input
	void SamplingPreProcessing ();

	/// creates the output mesh triangulation and initial vertices points from the clustering
	void BuildDelaunayTriangulation ();

	/// Adds polygons and vertices to fix boundaries
	void FixMeshBoundaries ();

	/// Experimental...
	void FixClusteringToVoronoi ();

	/// The window used to render the coarsened model
	RenderWindow *OutputMeshWindow;
	
	/// The window used to render the curvature indicator
	RenderWindow *IndicatorWindow;

	/// The output triangulation
	vtkSurface *Output;

	/// After the clustering is done, this method returns a list of Clusters adjacent to a specific dual Item
	/// Note that the Null cluster is not included in the list
	void GetDualItemNeighbourClusters (vtkIdType Item, vtkIdList * List);

	/// Adds a face in the coarsened triangulation
	vtkIdType AddFace (vtkIdType v1, vtkIdType v2, vtkIdType v3);

	/// Once the coarsened triangulation has been constructed, this method projects its vertices
	/// on the original surface, to make the approximation better
	void AdjustRemeshedGeometry ();

	/// Checks wether the input vtkSurface has to be subdivided or not before clustering 
	void CheckSubsamplingRatio ();

	/// Checks whether every output vertex is manifold
	/// returns the number of vertices with issues.
	int DetectNonManifoldOutputVertices ();

	/// the parameter storing the minimun subsampling ratio.
	/// if the actual subsampling ration is below, the input mesh will be subdivided accordingly
	/// default value is 10
	int SubsamplingThreshold;

	/// this pointer stores the input mesh if the mesh was subdivided before clustering
	vtkSurface *OriginalInput;

	/// this value stores the number of processed subdivisions
	int NumberOfSubdivisionsBeforeClustering;

	/// this array stores the parent-child informations (2 ints for each vertex: its two parents)
	/// it is used only when the mesh is bubdivided before simplification, to interpolate the 
	// curvature indicator;
	vtkIntArray *VerticesParent1;
	vtkIntArray *VerticesParent2;

	/// this option can be set to skip curvature computation and/or clustering by reading precomputed files
	int FileLoadSaveOption;

	// Experimental (Nevermind)
	int EdgeOptimization;

	// flag to enable addition of polygons and points to fix the mesh boundaries
	int BoundaryFixing;

	// name of custom imageData file giving user-defined density info
	char* InputDensityFile;

	double MaxCustomDensity;
	double MinCustomDensity;
	double CustomDensityMultiplicationFactor;
	bool ForceManifold;

};

template < class Metric >
int vtkDiscreteRemeshing < Metric >::DetectNonManifoldOutputVertices() {

	cout << "Starting detection of non-manifold vertices" << endl;

	// initialize items contained in clusters
	std::vector< vtkIdList *> ClusterItems;
	ClusterItems.resize( this->NumberOfClusters );

	for (vtkIdType i = 0; i != this->NumberOfClusters;i++ ) {

		ClusterItems[ i ] = 0;
		this->IsClusterFreezed->SetValue( i,1 );

	}

	// create a list of items contained in each cluster
	int numItems = this->GetNumberOfItems();
	int NumberOfUncorrectlyAssociatedItems = 0;

	for ( vtkIdType i = 0; i != numItems; i++ ) {

		int cluster = this->Clustering->GetValue( i );

		if( ( cluster < 0 ) || ( cluster >= this->NumberOfClusters ) )
			NumberOfUncorrectlyAssociatedItems++;
		else {

			if ( ClusterItems[ cluster ] == 0 )
				ClusterItems[ cluster ] = vtkIdList::New();

			ClusterItems[ cluster ]->InsertNextId(i);

		}

	}

	if ( NumberOfUncorrectlyAssociatedItems != 0 )
		cout << NumberOfUncorrectlyAssociatedItems << " uncorrectly associated items" << endl;

	int NumberOfTopologyIssues = 0;
	vtkIdList *ClustersWithIssues = vtkIdList::New();
	vtkIdList *CList = vtkIdList::New();

	for ( int cluster = 0; cluster != this->NumberOfClusters; cluster++ ) {

		if ( this->Output->IsVertexManifold( cluster ) ) continue;
		cout<<"Cluster "<<cluster<<" is non manifold"<<endl;
		vtkIdList *Items = ClusterItems[ cluster ];
		bool problem = true;

		if ( Items != 0) {
			cout << Items->GetNumberOfIds() << " items inside" << endl;
			if ( Items->GetNumberOfIds() == 1 )	{

				if ( this->Input->IsVertexManifold( Items->GetId( 0 ) ) != 1 ) {

					problem = false;
					cout << "discarding this topology issue as the input mesh also has a topology issue" << endl;

				}

			}

		} else {
			cout << ".... but empty. Skipping" << endl;
			continue;
		}

		if ( problem ) {

			NumberOfTopologyIssues++;
			ClustersWithIssues->InsertNextId(cluster);
/*			vtkIdList *CItems=ClusterItems[Cluster];
			int Size=0;
			if (CItems!=0)
				Size=CItems->GetNumberOfIds();

			cout<<"Vertex "<<Cluster<<" is non manifold. Cluster size : "<<Size<<endl;*/

			//unfreeze this cluster and its neighbours
			this->IsClusterFreezed->SetValue( cluster, 0 );
			this->Output->GetVertexNeighbours( cluster, CList );
			for ( int i = 0; i != CList->GetNumberOfIds(); i++ )
				this->IsClusterFreezed->SetValue( CList->GetId( i ), 0 );

		}

	}

	CList->Delete();
	vtkIdList *FList = vtkIdList::New();
	vtkIdList *IList = vtkIdList::New();

	for ( int i = 0; i != ClustersWithIssues->GetNumberOfIds(); i++ ) {

		vtkIdType cluster = ClustersWithIssues->GetId( i );
		vtkIdList *Items = ClusterItems[ cluster ];
		if ( Items == 0 ) {
			cout << "Warning : cluster " << cluster << " seems empty!" << endl;
			//continue;
		}
		int newSize = this->NumberOfClusters + 1;
		this->Clusters.resize( newSize );
		vtkIdType newCluster = this->NumberOfClusters;
		this->MetricContext.ResetCluster( &this->Clusters[ newCluster ] );
		ClusterItems.resize( newSize );
		ClusterItems[ newCluster ] = 0;
		this->IsClusterFreezed->Resize( newSize );
		this->IsClusterFreezed->SetValue( newCluster, false );
		this->ClustersLastModification.resize( newSize );
		this->ClustersLastModification[ newCluster ] = this->NumberOfLoops;
		this->ClustersSizes->Resize( newSize );
		this->ClustersSizes->SetValue( newCluster, 0 );
//		bool debug = newCluster == 22496 ? true : false;
//		if (debug ) cout<<"Adding cluster "<<newCluster<<" near cluster "<<Cluster<<endl;
		this->NumberOfClusters++;

		if ( Items->GetNumberOfIds() > 1 ) {

//			if (debug) cout << Items->GetNumberOfIds() << " items in cluster" << endl;
			vtkIdType ItemToMove = Items->GetId( 0 );
			this->Clustering->SetValue( ItemToMove,newCluster );
			ClusterItems[ newCluster ] = vtkIdList::New();
			ClusterItems[ newCluster ]->InsertNextId( ItemToMove );
			Items->DeleteId( ItemToMove );
			continue;

		}

//		if (debug) cout << "Only one item in cluster" << endl;
		// the cluster has only one item. Pick a neighbour item
		vtkIdType item = Items->GetId( 0 );
//			cout<<"Cluster "<<Cluster<<" contains only vertex "<<Item<<" of valence "
//			<<this->Input->GetValence(Item)<<endl;
		this->Input->GetVertexNeighbourFaces( item, FList );
/*			cout<<"Neighbour faces: "<<endl;
		for (int j=0;j!=FList->GetNumberOfIds();j++)
		{
			vtkIdType v1,v2,v3;
			this->Input->GetFaceVertices(FList->GetId(j),v1,v2,v3);
			cout<<FList->GetId(j)<<" : vertices "<<v1<<" "<<v2<<" "<<v3<<endl;
		}*/
		this->GetItemNeighbours( item, IList );
		bool found = false;

		for ( int j = 0; j < IList->GetNumberOfIds(); j++ ) {

			vtkIdType Neighbour = IList->GetId(j);
			vtkIdType NeighbourCluster = this->Clustering->GetValue( Neighbour );

			if ( ClusterItems[ NeighbourCluster ]->GetNumberOfIds() > 1 ) {

				this->Clustering->SetValue( Neighbour, newCluster );
				ClusterItems[ NeighbourCluster ]->DeleteId( Neighbour );
//					if ( NeighbourCluster == 22496 ) cout<<"Took 1 item from cluster "<<NeighbourCluster<<endl;
				ClusterItems[ newCluster ] = vtkIdList::New();
				ClusterItems[ newCluster ]->InsertNextId( Neighbour );
				found=true;
				break;
			}

		}

		if (!found) {

				cout<<"Could not find a place to add cluster "<<newCluster
				<<" near cluster "<<cluster<<endl;
			for ( int j = 0; j < IList->GetNumberOfIds(); j++ ) {

				vtkIdType Neighbour = IList->GetId( j );
				vtkIdType NeighbourCluster = this->Clustering->GetValue( Neighbour );

				if ( ClusterItems[ NeighbourCluster ]->GetNumberOfIds() > 1 ) {

					this->Clustering->SetValue( Neighbour, newCluster );
					ClusterItems[ NeighbourCluster ]->DeleteId( Neighbour );
					ClusterItems[ newCluster ]=vtkIdList::New();
					ClusterItems[ newCluster ]->InsertNextId( Neighbour );
					found = true;
					break;

				} else {
//						cout<<"Neighbour : "<<NeighbourCluster<<" has "
//							<<ClusterItems[NeighbourCluster]->GetNumberOfIds()<<" items"<<endl;
				}
			}

			this->Snapshot();

		}

		if ( !found ) {

			// no freeable nearby item found. Cancel new cluster creation
			newSize--;
			this->Clusters.resize( newSize );
			ClusterItems.resize( newSize );
			this->IsClusterFreezed->Resize( newSize );
			this->ClustersLastModification.resize( newSize );
			this->ClustersSizes->Resize( newSize );
			this->NumberOfClusters--;

		}

	}

	if ( this->Window ) 
		this->Window->DisplayRandomColors( this->NumberOfClusters + 1 );

	// free memory
	ClustersWithIssues->Delete();
	IList->Delete();
	for ( vtkIdType i = 0; i != this->NumberOfClusters; i++ )
		if ( ClusterItems[ i ] != 0 ) ClusterItems[ i ]->Delete();

	return NumberOfTopologyIssues;
}

template < class Metric >
	void vtkDiscreteRemeshing < Metric >::FixClusteringToVoronoi ()
{
/*
	int i,j;
	int Cluster,NearestCluster,Neighbour;
	double Distance,DistanceMin,ItemCoordinates[3],ClusterCoordinates[3];
	vtkIdList *NeighbourClusters=vtkIdList::New();

	for (i=0;i<this->GetNumberOfItems();i++)
	{
		Cluster=this->Clustering->GetValue(i);
		this->GetItemCoordinates(i,ItemCoordinates);
		this->GetClusterCentroid(Cluster,ClusterCoordinates);
		NearestCluster=Cluster;
		DistanceMin=vtkMath::Distance2BetweenPoints(ClusterCoordinates,ItemCoordinates);
		this->GetOutput()->GetVertexNeighbours(Cluster,NeighbourClusters);
		if (NeighbourClusters->GetNumberOfIds()>0)
		{
			cout<<NeighbourClusters->GetNumberOfIds()<<" ";
			for (j=0;j<NeighbourClusters->GetNumberOfIds();j++)
			{
				Neighbour=NeighbourClusters->GetId(j);
				this->GetClusterCentroid(Neighbour,ClusterBarycenter);
				Distance=vtkMath::Distance2BetweenPoints(ClusterBarycenter,TriangleBarycenter);
				if (Distance<DistanceMin)
				{
					DistanceMin=Distance;
					NearestCluster=Neighbour;
				}
			}
			this->Clustering->SetValue(i,NearestCluster);
		}
	}
	NeighbourClusters->Delete();
	this->ReComputeStatistics();
	this->Clustering->Modified();
	*/

}
template < class Metric >
	void vtkDiscreteRemeshing < Metric >::FixMeshBoundaries ()
{
	// until now, this method works only when clustering vertices
	if (this->ClusteringType == 0)
		return;

	if (this->BoundaryFixing == 0)
		return;

	vtkIdType i;
	vtkIdType *EdgesNewPoint = new vtkIdType[this->Input->GetNumberOfEdges ()];

	for (i = 0; i < this->GetInput ()->GetNumberOfEdges (); i++)
	{
		EdgesNewPoint[i] = -2;
	}


	vtkIdType f1, f2, v1, v2, c1, c2;
	double P[3], P1[3], P2[3];

	for (i = 0; i < this->GetInput ()->GetNumberOfEdges (); i++)
	{

		this->Input->GetEdgeVertices (i, v1, v2);
		this->Input->GetEdgeFaces (i, f1, f2);
		if ((f2 < 0) && (f1 >= 0))
		{
			// the edge is a boundary edge
			c1 = this->Clustering->GetValue (v1);
			c2 = this->Clustering->GetValue (v2);

			if (c1 != c2)
			{
				this->Input->GetPointCoordinates (v1, P1);
				this->Input->GetPointCoordinates (v2, P2);

				P[0] = 0.5 * (P1[0] + P2[0]);
				P[1] = 0.5 * (P1[1] + P2[1]);
				P[2] = 0.5 * (P1[2] + P2[2]);
				// create a new vertex in the coarsened model
				EdgesNewPoint[i]=this->Output->AddVertex (P);

				// create a new triangle (NewVertex,Cluster1,Cluster2);
				this->Output->AddFace (EdgesNewPoint[i],c1,c2);
			}
		}
	}

	vtkIdType Cluster;
	vtkIdType Edge;

	vtkIdList *EList = vtkIdList::New ();
	vtkIdType NumberOfEdges, *Edges;
	for (i = 0; i < this->GetInput ()->GetNumberOfEdges (); i++)
	{
		if (EdgesNewPoint[i] < 0)
		{
			// the edge was not visited previously 
			this->Input->GetEdgeFaces (i, f1, f2);
			if ((f2 < 0) && (f1 >= 0))
			{
				// the edge is a boundary edge 
				// with two vertices belonging to the same cluster

				this->Input->GetEdgeVertices (i, v1, v2);
				Cluster = this->Clustering->GetValue (v1);

				std::queue < vtkIdType >EQueue;
				EQueue.push (i);
				EList->Reset ();
				while (EQueue.size ())
				{
					Edge = EQueue.front ();
					EQueue.pop ();
					switch (EdgesNewPoint[Edge])
					{
					case -2:
						// this edge was never visited, so continue conquest
						EdgesNewPoint[Edge] = -1;
						this->Input->
							GetEdgeVertices (Edge,
									 v1,
									 v2);
						this->Input->
							GetVertexNeighbourEdges
							(v1, NumberOfEdges,
							 Edges);
						vtkIdType j;
						for (j = 0; j < NumberOfEdges;
						     j++)
						{
							this->Input->
								GetEdgeFaces
								(Edges[j], f1,
								 f2);
							if (f2 < 0)
								EQueue.push
									(Edges
									 [j]);
						}
						this->Input->
							GetVertexNeighbourEdges
							(v2, NumberOfEdges,
							 Edges);
						for (j = 0; j < NumberOfEdges;
						     j++)
						{
							this->Input->
								GetEdgeFaces
								(Edges[j], f1,
								 f2);
							if (f2 < 0)
								EQueue.push
									(Edges
									 [j]);
						}
					case -1:
						break;
					default:
						EList->InsertUniqueId (Edge);
					}
				}

				if (EList->GetNumberOfIds () >= 2)
				{
					this->Output->AddFace (Cluster,
							       EdgesNewPoint
							       [EList->
								GetId (0)],
							       EdgesNewPoint
							       [EList->
								GetId (1)]);
				}
			}
		}
	}

	// fix when two neighbour edges actually have several clusters
	for (i = 0; i < this->Input->GetNumberOfPoints (); i++)
	{
		if (this->Input->GetNumberOfBoundaries (i) > 0)
		{
			EList->Reset ();
			this->Input->GetVertexNeighbourEdges (i, NumberOfEdges, Edges);
			vtkIdType j;
			for (j = 0; j < NumberOfEdges; j++)
			{
				Edge = Edges[j];
				if (EdgesNewPoint[Edge] >= 0)
				{
					EList->InsertUniqueId (EdgesNewPoint[Edge]);
				}
			}
			if (EList->GetNumberOfIds () > 1)
			{
				this->Output->AddFace (this->Clustering->GetValue (i),EList->GetId (0),EList->GetId (1));
			}
		}
	}

	EList->Delete ();
	delete[]EdgesNewPoint;
}


template < class Metric >
void vtkDiscreteRemeshing < Metric >::SamplingPreProcessing () {

	int i;
	int Compute = 0;
	vtkDoubleArray *CellsIndicators;
	vtkDoubleArray *CustomIndicatorColors = vtkDoubleArray::New ();
	CustomIndicatorColors->SetNumberOfValues (this->GetNumberOfItems ());

	if (this->MetricContext.IsCurvatureIndicatorNeeded () == 1)
	{
		if (this->FileLoadSaveOption)
		{
			cout << "Do you want to compute the curvature	indicator? (0:No 1:Yes)	";
			cin >> Compute;
		}

		vtkPolyData *PrincipalDirectionsPolyData = 0;
		vtkDataArrayCollection *CurvatureCollection=0;
		if ((Compute == 1) || (this->FileLoadSaveOption == 0))
		{
			if (this->InputDensityFile)
			{
				vtkImageReader2 *Reader;
				if(strstr (this->InputDensityFile,".mnc") != NULL)
				{
					Reader=vtkMINCImageReader::New();
					((vtkMINCImageReader*) Reader)->RescaleRealValuesOn();
				}
				else
					Reader=vtkMetaImageReader::New();

				Reader->SetFileName(this->InputDensityFile);
				Reader->Update();
				vtkImageData *Density=Reader->GetOutput();
				CellsIndicators=vtkDoubleArray::New();
				int NumberOfItems=this->GetNumberOfItems();
				CellsIndicators->SetNumberOfValues(NumberOfItems);
				CurvatureCollection=vtkDataArrayCollection::New();
				CurvatureCollection->AddItem(CellsIndicators);

				for (int i=0;i<NumberOfItems;i++)
				{
					double Coord[3];
					int ijk[3];
					double pcoords[3];
					this->GetItemCoordinates(i,Coord);
					Density->ComputeStructuredCoordinates(Coord,ijk,pcoords);
					double Value=Density->GetScalarComponentAsDouble(ijk[0],ijk[1],ijk[2],0);
					Value=this->MaxCustomDensity-Value*this->CustomDensityMultiplicationFactor;
					if (Value<this->MinCustomDensity)
						Value=this->MinCustomDensity;
					CellsIndicators->SetValue(i,Value);
				}
				Reader->Delete();
			}
			else
			{
				vtkCurvatureMeasure *Curvature =vtkCurvatureMeasure::New ();
				if (this->OriginalInput)
					Curvature->SetInputData (this->OriginalInput);
				else
					Curvature->SetInputData (this->Input);
			
				Curvature->SetComputationMethod (1);
				Curvature->SetElementsType (this->ClusteringType);
				Curvature->SetComputePrincipalDirections (this->MetricContext.IsPrincipalDirectionsNeeded());

				CurvatureCollection=Curvature->GetCurvatureIndicator();
				CurvatureCollection->Register(this);
			
				CellsIndicators =(vtkDoubleArray *) CurvatureCollection->GetItem (0);

				if (this->FileLoadSaveOption != 0)
				{
					fstream CurvatureOutput;
					CurvatureOutput.open ("curvature.dat", std::ofstream::out | std::ofstream::trunc | std::ios::binary);
					int NumberOfItems;
					vtkSurface *InputSurface;
					if (this->OriginalInput)
						InputSurface=this->OriginalInput;
					else
						InputSurface=this->Input;
				
					if (this->ClusteringType == 0)
						NumberOfItems =InputSurface->GetNumberOfCells ();
					else
						NumberOfItems =InputSurface->GetNumberOfPoints ();

					for (i = 0; i < NumberOfItems; i++)
					{
						double value;
						value = CellsIndicators->GetValue (i);
						CurvatureOutput.write ((char *) &value,sizeof (double));
					}
					CurvatureOutput.close ();
				}

				if ((this->MetricContext.IsPrincipalDirectionsNeeded () == 1)&&(this->Display!=0))
					PrincipalDirectionsPolyData =Curvature->GetPrincipalDirectionsPolyData ();

				Curvature->Delete ();
			}
		}
		else
		{
			CellsIndicators = vtkDoubleArray::New ();
			CellsIndicators->SetNumberOfValues (this->GetNumberOfItems());

			fstream CurvatureInput;
			CurvatureInput.open ("curvature.dat",std::ofstream::in | std::ios::binary);

			for (i = 0; i < this->GetNumberOfItems (); i++)
			{
				double value;
				CurvatureInput.read ((char *) &value,sizeof (double));
				CellsIndicators->SetValue (i, value);
			}
			CurvatureInput.close ();
		}

		// now we have to interpolate the curvature measure when the input mesh was subdivided
		if (this->NumberOfSubdivisionsBeforeClustering != 0)
		{
			cout << "Interpolating...";
			vtkDoubleArray *CellsIndicators2 =vtkDoubleArray::New ();
			CellsIndicators2->SetNumberOfValues (this->GetNumberOfItems());

			if (this->ClusteringType == 0)
			{
				// Interpolating faces indicators
				int NumberOfChildrenFacesPerFace;
				NumberOfChildrenFacesPerFace =4 << (2 *(this->NumberOfSubdivisionsBeforeClustering- 1));
				int CellIndex = 0;
				double CurvatureValue;
				double Ratio =
					(double) NumberOfChildrenFacesPerFace;
				int j;
				for (i = 0;i < OriginalInput->GetNumberOfCells ();i++)
				{
					CurvatureValue =CellsIndicators->GetValue (i) / Ratio;
					for (j = 0;j < NumberOfChildrenFacesPerFace;j++)
					{
						CellsIndicators2->SetValue (CellIndex,CurvatureValue);
						CellIndex++;
					}
				}
			}
			else
			{
				// Interpolating vertices indicators
				for (i = 0; i < OriginalInput->GetNumberOfPoints ();i++)
				{
					CellsIndicators2->SetValue (i,CellsIndicators->GetValue(i));

				}
				for (i = OriginalInput->GetNumberOfPoints ();i < this->Input->GetNumberOfPoints ();i++)
				{
					CellsIndicators2->SetValue (i, 
						(						
						CellsIndicators2->GetValue(this->VerticesParent1->GetValue(i)) +
						CellsIndicators2->GetValue(this->VerticesParent2->GetValue(i))
						) *0.5);
				}
			}
			CurvatureCollection->ReplaceItem(0,CellsIndicators2);
			CellsIndicators2->Delete();
//			CellsIndicators->Delete();
			CellsIndicators = CellsIndicators2;

			cout << ".... Done" << endl;
		}
		double Range[2];
		CellsIndicators->GetRange (Range);
		cout<<"Indicators Range : "<<Range[0]<<"  "<<Range[1]<<endl;
		double MaxIndicatorColor = 0;
		double MeanIndicator = 0;

		// Compute Average and Maximum indicator
		for (i = 0; i < this->GetNumberOfItems (); i++)
		{
			CustomIndicatorColors->SetValue (i,pow(CellsIndicators->GetValue (i),
						this->MetricContext.GetGradation ()));
			if (MaxIndicatorColor <CustomIndicatorColors->GetValue (i))
				MaxIndicatorColor =CustomIndicatorColors->GetValue (i);
			MeanIndicator += CustomIndicatorColors->GetValue (i);
		}

		MeanIndicator /= (double) this->GetNumberOfItems ();

		double ClampingFactor = 100.0;
		if (MaxIndicatorColor > ClampingFactor * MeanIndicator)
			MaxIndicatorColor = ClampingFactor * MeanIndicator;
		for (i = 0; i < this->GetNumberOfItems (); i++)
		{
			if (CustomIndicatorColors->GetValue (i) >MaxIndicatorColor)
				CustomIndicatorColors->SetValue (i,MaxIndicatorColor);
			CustomIndicatorColors->SetValue (i,CustomIndicatorColors->GetValue (i)/MaxIndicatorColor);

		}

		this->MetricContext.SetCurvatureInfo (CurvatureCollection);
		CurvatureCollection->Delete();

		if (this->Display)
		{
			RenderWindow *Window;
			vtkPolyData *Mesh2 = vtkPolyData::New ();
			Window = RenderWindow::New ();

			Mesh2->ShallowCopy (this->Input);

			if (this->ClusteringType == 0)
			{
				Mesh2->GetPointData()->SetScalars (0);
				Mesh2->GetCellData()->SetScalars (CustomIndicatorColors);
			}
			else
			{
				Mesh2->GetPointData ()->SetScalars (CustomIndicatorColors);
				Mesh2->GetCellData ()->SetScalars (0);
			}

			Window->SetInputData (Mesh2);
			Mesh2->Delete();
			if (PrincipalDirectionsPolyData)
			{
				Window->SetInputEdges (PrincipalDirectionsPolyData);
				PrincipalDirectionsPolyData->Delete();
			}	

			vtkLookupTable *bwLut = vtkLookupTable::New ();
			bwLut->SetTableRange (0, 1);

			bwLut->SetSaturationRange (0, 0);	// Black     and     white
			bwLut->SetHueRange (0, 0);
			bwLut->SetValueRange (0, 1);
			bwLut->SetScaleToLog10 ();
			bwLut->Build ();
			Window->SetLookupTable (bwLut);

			Window->Render ();
			Window->SetWindowName ("Curvature Indicator");
			if (this->AnchorRenderWindow)
				Window->AttachToRenderWindow (this->AnchorRenderWindow);
			else
				this->AnchorRenderWindow = Window;
			Window->Interact ();
			if (this->IndicatorWindow)
				this->IndicatorWindow->Delete();
			this->IndicatorWindow=Window;

			bwLut->Delete ();
		}
	}
	CustomIndicatorColors->Delete();
}

template < class Metric >
void vtkDiscreteRemeshing < Metric >::CheckSubsamplingRatio () {

	vtkSurface *Levels[ 100 ];
	Levels[ 0 ] = this->Input;

	if ( this->ClusteringType == 1 ) {

		if ( !this->VerticesParent1 )
			this->VerticesParent1 = vtkIntArray::New ();
		if (!this->VerticesParent2)
			this->VerticesParent2 = vtkIntArray::New ();

	}

	while ( Levels[ NumberOfSubdivisionsBeforeClustering ]->GetNumberOfPoints()
		< this->SubsamplingThreshold * this->NumberOfClusters ) {

		if ( this->ConsoleOutput ) cout << "Subdividing mesh" << endl;
		Levels[ NumberOfSubdivisionsBeforeClustering + 1 ] =
			Levels[ NumberOfSubdivisionsBeforeClustering ]->
			Subdivide( this->VerticesParent1, this->VerticesParent2 );
		NumberOfSubdivisionsBeforeClustering++;

	}

	if ( NumberOfSubdivisionsBeforeClustering ) {

		this->OriginalInput = this->Input;
		this->Input = Levels[ NumberOfSubdivisionsBeforeClustering ];
		for ( int i = 1; i < NumberOfSubdivisionsBeforeClustering; i++ )
			Levels[i]->Delete ();

	}

}

template < class Metric > void vtkDiscreteRemeshing < Metric >::Remesh ()
{
	int Compute = 1;
	this->CheckSubsamplingRatio ();
	this->SamplingPreProcessing ();

	if ( this->ConsoleOutput )
		cout << "Input mesh: " << this->GetInput ()->GetNumberOfPoints()
			<< " vertices	and	" << this->GetInput ()->GetNumberOfCells()
			<< " faces" << endl;

	if ( this->FileLoadSaveOption ) {

		cout << "Do you want to compute the clustering? (0:NO	1:Yes) ";
		cin >> Compute;

	}

	if ( ( Compute == 1) || ( this->FileLoadSaveOption == 0 ) ) {

		this->ProcessClustering ();

		if ( this->FileLoadSaveOption == 1 ) {

			fstream ClusteringOutput;
			ClusteringOutput.open( "clustering.dat",
				std::ofstream::out | std::ofstream::trunc | std::ios::binary );

			for ( int i = 0; i < this->GetNumberOfItems (); i++ ) {

				int value = this->Clustering->GetValue( i );
				ClusteringOutput.write ( (char *) &value, sizeof (int) );
			}

			ClusteringOutput.close ();

		}

	} else {

		this->Init ();
		fstream ClusteringInput;
		ClusteringInput.open( "clustering.dat", std::ofstream::in | std::ios::binary );

		for ( int i = 0; i < this->GetNumberOfItems (); i++) {

			int value;
			ClusteringInput.read( (char *) &value, sizeof ( int ) );
			this->Clustering->SetValue( i, value );

		}

		ClusteringInput.close();
		this->ReComputeStatistics();

	}

	this->BuildDelaunayTriangulation ();
	if ( this->ConsoleOutput ) this->Output->DisplayMeshProperties();

	if ( this->ForceManifold ) {

		while ( int NumberOfIssues = this->DetectNonManifoldOutputVertices() ) {

			cout << NumberOfIssues << " topology issues, restarting minimization" << endl;
			this->ConnexityConstraint = 0;
			this->MinimizeEnergy();
			this->BuildDelaunayTriangulation();

		}

		if ( this->ConsoleOutput ) this->Output->DisplayMeshProperties ();

	}

}
template < class Metric >
void vtkDiscreteRemeshing <Metric >::GetDualItemNeighbourClusters ( vtkIdType item,vtkIdList * List ) {

	if (this->ClusteringType == 0) {

		vtkIdType cluster;
		List->Reset ();
		vtkIdList *VerticesList = vtkIdList::New ();
		this->GetInput ()->GetVertexNeighbourFaces( item,VerticesList );

		for ( vtkIdType i = 0; i < VerticesList->GetNumberOfIds (); i++ ) {

			cluster = this->Clustering->GetValue( VerticesList->GetId( i ) );
			if ( cluster != this->NumberOfClusters ) List->InsertUniqueId( cluster );

		}

		VerticesList->Delete ();

	} else {

		vtkIdType *Vertices;
		vtkIdType NumberOfVertices;
		List->Reset ();
		this->Input->GetFaceVertices( item, NumberOfVertices, Vertices );

		for( vtkIdType i = 0; i < NumberOfVertices; i++ ) {

			vtkIdType cluster = this->Clustering->GetValue( Vertices[ i ] );
			if ( cluster < this->NumberOfClusters )	List->InsertUniqueId( cluster );

		}

	}

}

template < class Metric >
vtkIdType vtkDiscreteRemeshing < Metric >::AddFace( vtkIdType v1, vtkIdType v2, vtkIdType v3) {

	if ( ( v1 == v2 ) || ( v1 == v3 ) || ( v2 == v3 ) ) return (-1);
	if ( this->Output->IsFace( v1, v2, v3 ) < 0 )
		return ( this->Output->AddFace( v1, v2, v3 ) );
	else return -1;

}

template < class Metric >
void vtkDiscreteRemeshing < Metric >::BuildDelaunayTriangulation () {

	double P[3];
	vtkIdList *CList = vtkIdList::New ();
	vtkIdType v1, v2, v3, f1, f2;
	if ( this->Output != 0 ) this->Output->Delete();
	this->Output = vtkSurface::New();

	// Find the first non-empty cluster and put its Id in Valid
	int Valid = 0;

	for ( int i = 0; i < this->GetNumberOfClusters (); i++ )	{

		if ( this->ClustersSizes->GetValue( i ) > 0 ) {

			Valid = i;
			break;

		}

	}

	// We compute the vertices as inertia centers of each Cluster
	for ( int i = 0; i < this->NumberOfClusters; i++ ) {

		if ( this->ClustersSizes->GetValue (i) == 0 )
			this->MetricContext.GetClusterCentroid ( &this->Clusters[ Valid ], P );
		else
			this->MetricContext.GetClusterCentroid ( &this->Clusters[ i ], P );

		this->Output->AddVertex ( P[ 0 ], P[ 1 ], P[ 2 ] );

	}

	for ( int i = 0; i < this->GetNumberOfDualItems (); i++ ) {

		this->GetDualItemNeighbourClusters ( i, CList );

		if ( this->ClusteringType == 0 ) {

			if ( CList->GetNumberOfIds () > 3 ) {

				if ( ( this->GetInput ()->GetNumberOfBoundaries ( i ) == 0)
				    && ( CList->GetNumberOfIds () < 600 ) ) {

					CList->Reset ();
					vtkIdType e1 = this->GetInput ()->GetFirstEdge( i );
					this->GetInput ()->GetEdgeVertices( e1, v1, v2 );
					if ( v1 == i ) v1 = v2;
					vtkIdType v_init = v1;
					this->GetInput ()->GetEdgeFaces( e1, f1, f2 );
					v2 = this->GetInput ()->GetThirdPoint( f1, i, v1 );
					v1 = -1;
					vtkIdType type_last = -1;
					int Valence = this->GetInput()->GetValence( i );

					while ( Valence >= 0 ) {

						Valence--;
						if ( v1 == v_init )	break;
						vtkIdType type = this->Clustering->GetValue( f1 );

						if ( type != type_last ) {

							CList->InsertNextId( type );
							type_last = type;

						}

						this->GetInput()->Conquer( f1, i, v2, f2, v3 );
						if ( f2 < 0 ) break;
						f1 = f2;
						v1 = v2;
						v2 = v3;

					}

				}

			}

			vtkIdType n = CList->GetNumberOfIds ();

			if ( n >= 3 ) {

				v1 = CList->GetId( 0 );
				if ( v1 == CList->GetId( n - 1 ) ) n--;

				for ( int j = 0; j < n - 2; j++)
					this->AddFace( v1, CList->GetId( j + 1 ), CList->GetId( j + 2 ) );

			}

		} else {

			if ( CList->GetNumberOfIds() == 3 )
				this->AddFace( CList->GetId( 0 ), CList->GetId( 1 ), CList->GetId( 2 ) );

			if ( CList->GetNumberOfIds () == 4 ) {

				this->AddFace( CList->GetId( 0 ), CList->GetId( 1 ), CList->GetId( 2 ) );
				this->AddFace( CList->GetId( 0 ), CList->GetId( 2 ), CList->GetId( 3 ) );

			}

		}

	}

	CList->Delete();
	this->FixMeshBoundaries();

	if ( this->ForceManifold ) {

		// add non manifold edges
		for ( vtkIdType i = 0; i < this->GetNumberOfEdges(); i++ ) {

			vtkIdType I1,I2;
			vtkIdType C1,C2;
			this->GetEdgeItems( i, I1, I2 );
			C1 = this->Clustering->GetValue( I1 );
			C2 = this->Clustering->GetValue( I2 );

			if ( ( C1 != C2 )
				&& ( C1 >= 0 ) && ( C1 < this->NumberOfClusters )
				&& ( C2 >= 0 ) && ( C2 < this->NumberOfClusters ) ) {

				if ( this->Output->IsEdge( C1, C2 ) < 0 )
					this->Output->AddEdge(C1,C2);
			}
		}
	}


	if ( this->Display ) {

		if (this->OutputMeshWindow == 0 )
			this->OutputMeshWindow = RenderWindow::New ();

		this->OutputMeshWindow->SetInputData( this->Output );
		this->OutputMeshWindow->DisplayInputEdges ();
		this->OutputMeshWindow->Render();
		this->OutputMeshWindow->SetWindowName ("Coarsened model");
		this->OutputMeshWindow->SetLookupTable();
		if ( this->AnchorRenderWindow )
			this->OutputMeshWindow->AttachToRenderWindow( this->AnchorRenderWindow );
		this->OutputMeshWindow->Interact();

	}

}

template < class Metric >
void vtkDiscreteRemeshing < Metric >::AdjustRemeshedGeometry() {

	int i, Cluster;
	vtkMath *Math = vtkMath::New ();

	double *Distances = new double[this->GetNumberOfClusters ()];
	double P[3], Pm[3], distance;

	vtkPoints *Points = vtkPoints::New ();
	Points->DeepCopy (this->Output->GetPoints ());

	for (i = 0; i < this->GetNumberOfClusters (); i++)
		Distances[i] = 100000000;

	for (i = 0; i < this->GetNumberOfItems (); i++)
	{
		Cluster = this->Clustering->GetValue (i);
		this->GetItemCoordinates (i, P);

		Points->GetPoint (Cluster, Pm);
		distance = sqrt (Math->Distance2BetweenPoints (P, Pm));
		if (distance < Distances[Cluster])
		{
			Distances[Cluster] = distance;
			this->Output->SetPointCoordinates (Cluster, P);
		}
	}

	Math->Delete ();
	delete[] Distances;
	Points->Delete ();

}

template < class Metric >
vtkDiscreteRemeshing < Metric >::vtkDiscreteRemeshing () {

	this->BoundaryFixing = 0;
	this->EdgeOptimization = 0;
	this->AnchorRenderWindow = 0;
	this->FileLoadSaveOption = 0;
	this->OriginalInput = 0;
	this->VerticesParent1 = 0;
	this->VerticesParent2 = 0;
	this->SubsamplingThreshold = 10;
	this->NumberOfSubdivisionsBeforeClustering = 0;
	this->OutputMeshWindow = 0;
	this->IndicatorWindow = 0;
	this->InputDensityFile  =0;
	this->MaxCustomDensity = 1;
	this->MinCustomDensity = 0.1;
	this->CustomDensityMultiplicationFactor = 0.001;
	this->Output = 0;
	this->ForceManifold = false;

}

template < class Metric >
	vtkDiscreteRemeshing < Metric >::~vtkDiscreteRemeshing ()
{
	if (this->OriginalInput)
		this->OriginalInput->Delete();
		
	if (this->Output)
		this->Output->Delete();
	
	if (this->VerticesParent1)
		this->VerticesParent1->Delete();
		
	if (this->VerticesParent2)
		this->VerticesParent2->Delete();
	
	if (this->OutputMeshWindow)
		this->OutputMeshWindow->Delete();
	
	if (this->IndicatorWindow)
		this->IndicatorWindow->Delete();
}
#endif
