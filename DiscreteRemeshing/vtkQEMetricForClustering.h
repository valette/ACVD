/*=========================================================================

Program:   Quadric Error Metric for Clustering 
Module:    vtkQEMetricForClustering.h
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

#ifndef _VTKQEMETRICFORCLUSTERING_H_
#define _VTKQEMETRICFORCLUSTERING_H_


/// the Quadrics Error metrics defined by the paper: "Surface simplification using Quadric Error Metrics",
// Garland and Heckbert, Siggraph 97

#include <vtkTriangle.h>
#include <vtkMath.h>
#include <vtkDataArrayCollection.h>

#include "vtkQuadricTools.h"

class  vtkQEMetricForClustering {

public:

	void SetQuadricsOptimizationLevel( int L ) {

		this->QuadricsOptimizationLevel = L;

	}

	void SetConstrainedClustering( int C ) {

		cout << "Setting Metric Constraint to "<< C << endl;
		this->ActiveConstraints = C;

	}

	int IsCurvatureIndicatorNeeded() {

		return this->Gradation > 0 ? 1 : 0;

	}
	
	int IsPrincipalDirectionsNeeded() {

		return 0;

	}
	
	void SetCurvatureInfo( vtkDataArrayCollection *Info ) {

		if ( this->CustomWeights )
			this->CustomWeights->Delete();

		this->CustomWeights=( vtkDoubleArray* ) Info->GetItem( 0 );
		this->CustomWeights->Register( this->Object );

	}


	void SetGradation( double Gradation ) {

		this->Gradation = Gradation;

	}

	double GetGradation() { return this->Gradation; }

	// each item is represented by a quadric matrix
	struct Item {

		// 4x4 quadric in compact form, without the 10th useless coefficient
		double Quadric[ 9 ];

		// barycenter (multiplied by its weight)
		double Value[ 3 ];

		// weight
		double Weight;
	};

	struct Cluster {

		// Value of the energy term for this cluster (cached to increase speed by about 33%)
		double EnergyValue;

		// accumulated quadrics, without the 10th useless coefficient
		double SQuadric[ 9 ];

		// cluster best approximating point ( centroid, quadric-based or fixed vertex)
		double Centroid[ 3 ];

		// accumulated centroids (weighted)
		double SValue[ 3 ];

		// accumulated weights
		double SWeight;

		/// gives the number of quadric eigenvalues equal to zero (could be 0, 1 or 2)
		char RankDeficiency;

		Cluster () { AnchorItem = -1; }
		vtkIdType AnchorItem;
	};

	int GetClusterRankDeficiency( Cluster *C ) {

		return C->RankDeficiency;

	}

	void ResetItem( Item *I ) {

		for ( int i = 0; i < 9; i++ ) I->Quadric[ i ] = 0;
		for ( int i = 0; i < 3; i++ ) I->Value[ i ] = 0;
		I->Weight = 0;

	}

	double GetItemWeight( vtkIdType ItemId ) {

		return this->Items[ ItemId ].Weight;

	}

	void ComputeVertexQuadric( int Vertex ) {

		Item *I = this->Items + Vertex;
		double Weight = I->Weight;
		this->ResetItem( I );
		vtkIdList *FList=vtkIdList::New();
		Mesh->GetVertexNeighbourFaces(Vertex,FList);

		for ( int i = 0; i < FList->GetNumberOfIds(); i++ )
			vtkQuadricTools::AddTriangleQuadric( I->Quadric, Mesh, FList->GetId( i ), false );

		Mesh->GetPoint( Vertex, I->Value );
		for ( int i = 0; i < 3; i++) I->Value[ i ] = I->Value[ i ] * Weight;
		I->Weight = Weight;
		FList->Delete();

	}

	void ComputeTriangleQuadric( int Face ) {

		double x1[ 3 ], x2[ 3 ], x3[ 3 ];
		vtkIdType V1, V2, V3;
		Item *I = this->Items + Face;
		double Weight = I->Weight;
		this->ResetItem( I );
		vtkQuadricTools::AddTriangleQuadric( I->Quadric, Mesh, Face, false);
		Mesh->GetFaceVertices( Face, V1, V2, V3 );
		Mesh->GetPointCoordinates( V1, x1 );
		Mesh->GetPointCoordinates( V2, x2 );
		Mesh->GetPointCoordinates( V3, x3 );
		I->Weight = Weight;

		for ( int i = 0; i < 3; i++ )
			I->Value[ i ] = Weight * ( x1[ i ] + x2[ i ] + x3[ i ] ) / 3.0;

	}

	double GetClusterEnergy( Cluster *C ) {

		return C->EnergyValue;

	}

	void ComputeClusterEnergy( Cluster *C ) {

		if ( C->AnchorItem == -2 ) {

			C->EnergyValue = 1e100;
			return;

		}

		C->EnergyValue= (
			C->Centroid[ 0 ] * C->Centroid[ 0 ]
			+ C->Centroid[ 1 ]*C->Centroid[ 1 ]
			+ C->Centroid[ 2 ]*C->Centroid[ 2 ]) * C->SWeight
			-2.0 * vtkMath::Dot( C->Centroid, C->SValue );
	}

	void DeepCopy( Cluster *Source, Cluster *Destination ) {

		Destination->SWeight = Source->SWeight;

		for ( int i = 0; i < 3; i++ ) {

			Destination->SValue[ i ] = Source->SValue[ i ];
			Destination->Centroid[ i ] = Source->Centroid[ i ];

		}

		for ( int i = 0; i < 9; i++ )
			Destination->SQuadric[ i ] = Source->SQuadric[ i ];

		Destination->EnergyValue = Source->EnergyValue;
		Destination->RankDeficiency = Source->RankDeficiency;
		Destination->AnchorItem = Source->AnchorItem;

	}

	void Add( Cluster *Source, vtkIdType ItemId, Cluster *Destination ) {

		this->DeepCopy( Source, Destination );
		this->AddItemToCluster( ItemId, Destination );

	}

	void AddItemToCluster( vtkIdType ItemId, Cluster *C ) {

		Item *I = &this->Items[ ItemId ];
		for ( int i = 0; i < 3; i++ ) C->SValue[ i ] += I->Value[ i ];
		for ( int i = 0; i < 9; i++) C->SQuadric[ i ] += I->Quadric[ i ];
		C->SWeight+=I->Weight;
	}

	void Sub( Cluster *Source, vtkIdType ItemId, Cluster *Destination ) {

		this->DeepCopy( Source,Destination );
		this->SubstractItemFromCluster( ItemId, Destination );

	}

	void SubstractItemFromCluster( vtkIdType ItemId, Cluster *C) {

		Item *I = &this->Items[ ItemId ];
		for ( int i = 0; i < 3; i++ ) C->SValue[ i ] -= I->Value[ i ];
		for ( int i = 0; i < 9; i++ ) C->SQuadric[ i ] -= I->Quadric[ i ];
		C->SWeight -= I->Weight;

		if ( C->AnchorItem == ItemId ) C->AnchorItem = -2; // broken constraint

	}

	void GetClusterCentroid( Cluster *C, double *P ) {

		for ( int i = 0; i < 3; i++ ) P[ i ]=C->Centroid[ i ];
	}

	void ComputeClusterCentroid( Cluster *C ) {

		if ( C->AnchorItem >= 0 ) {

			Mesh->GetPointCoordinates( C->AnchorItem, C->Centroid );
			return;

		}

		for ( int i = 0; i < 3; i++) C->Centroid[ i ] = C->SValue[ i ] / C->SWeight;
		
		if ( !this->ActiveConstraints || !this->QuadricsOptimizationLevel )
			return;

		C->RankDeficiency = vtkQuadricTools::ComputeRepresentativePoint(
			C->SQuadric, C->Centroid, this->QuadricsOptimizationLevel );

	}

	void ResetCluster( Cluster *C ) {

		for ( int i = 0; i < 9; i++ ) C->SQuadric[ i ] = 0;
		for ( int i = 0; i < 3; i++ ) C->SValue[ i ]=0;
		for ( int i = 0; i < 3; i++ ) C->Centroid[ i ]=0;
		C->SWeight = 0;
		C->EnergyValue = 0;

	}

	void BuildMetric( vtkSurface *Mesh, int ClusteringType ) {

		this->Mesh = Mesh;

		if ( ClusteringType == 0 ) {
			// Items are triangles
			vtkIdType v1,v2,v3;
			double P1[ 3 ],P2[ 3 ],P3[ 3 ];
			Items = new Item[ Mesh->GetNumberOfCells() ];

			for ( vtkIdType i=0;i < Mesh->GetNumberOfCells(); i++ ) {

				Mesh->GetFaceVertices( i, v1, v2, v3 );
				Mesh->GetPoint( v1, P1 );
				Mesh->GetPoint( v2, P2 );
				Mesh->GetPoint( v3, P3 );

				Items[ i ].Weight = vtkTriangle::TriangleArea( P1, P2, P3 );
				if ( Gradation > 0 ) Items[ i ].Weight *=
					pow( CustomWeights->GetValue( i ), Gradation );

			}

			this->ClampWeights( Items, Mesh->GetNumberOfCells(), 10000 );
			for ( vtkIdType i = 0; i < Mesh->GetNumberOfCells(); i++ )
				this->ComputeTriangleQuadric( i );

		} else {

			// Items are vertices
			Items = new Item[ Mesh->GetNumberOfPoints() ];
			vtkDoubleArray *VerticesAreas = Mesh->GetVerticesAreas();
			for ( vtkIdType i = 0; i < Mesh->GetNumberOfPoints(); i++ ) {

				Items[i].Weight = VerticesAreas->GetValue(i);

				if ( Gradation > 0 ) Items[i].Weight *=
					pow( CustomWeights->GetValue( i ), Gradation );

			}

			this->ClampWeights( Items,Mesh->GetNumberOfPoints(), 10000 );
			for ( vtkIdType i = 0; i < Mesh->GetNumberOfPoints();i++ )
				this->ComputeVertexQuadric( i );

		}

	}

	// this method clamps the weights between AverageValue/Ratio and AverageValue*Ratio
	void ClampWeights( Item *Items, int NumberOfValues, double Ratio ) {

		double Average = 0;
		for ( int i = 0; i < NumberOfValues; i++ ) Average += Items[ i ].Weight;
		Average /= (double) NumberOfValues;
		double Min = Average/Ratio;
		double Max = Average*Ratio;

		for ( int i = 0; i < NumberOfValues; i++ ) {

			if ( Items[ i ].Weight > Max ) Items[ i ].Weight = Max;
			if ( Items[ i ].Weight < Min ) Items[ i ].Weight = Min;

		}

	}

	vtkQEMetricForClustering() {

		this->CustomWeights = 0;
		this->Gradation = 0;
		this->ActiveConstraints = 1;
		this->Object = vtkObject::New();
		this->QuadricsOptimizationLevel = 3;

	};

	~vtkQEMetricForClustering() {

		this->Object->Delete();
		if ( this->CustomWeights ) this->CustomWeights->Delete();

	};

private:

	vtkSurface *Mesh;
	int	ActiveConstraints;
	vtkDoubleArray *CustomWeights;
	int QuadricsOptimizationLevel;
	vtkObject *Object; // Dummy object used for registering the curvature indicators
	double Gradation;
	Item *Items; // array storing items

};

#endif
