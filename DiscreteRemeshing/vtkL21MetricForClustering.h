/*=========================================================================

Program:   L2,1 Metric for Clustering 
Module:    vtkL21MetricForClustering.h
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

#ifndef _VTKL21METRICFORCLUSTERING_H_
#define _VTKL21METRICFORCLUSTERING_H_

#include <vtkTriangle.h>
#include <vtkMath.h>
#include <vtkIdList.h>
#include "vtkSurface.h"

// This metric is defined in the paper "Variationnal shape Approximation", Cohen-Steiner et al., Siggraph 2004

class  vtkL21MetricForClustering 
{
public:

	int IsCurvatureIndicatorNeeded()
	{
		return (0);
	}

	void SetCustomWeights(vtkDoubleArray *Weights)
	{
		this->CustomWeights=Weights;
	}

	void SetGradation(double Gradation)
	{
		this->Gradation=Gradation;
	}

	double GetGradation()
	{
		return (this->Gradation);
	}

	// each item is represented by a quadric matrix
	struct Item
	{
		// the item normal vector, wheighted by its area
		double Normal[3];

		// the barycenter (multiplied by its weight)
		double Value[3];

		// the weight
		double Weight;
	};

	void ResetItem(Item *I) 
	{
		int i;
		for (i=0;i<3;i++)
		{
			I->Normal[i]=0;
			I->Value[i]=0;
		}
		I->Weight=0;
	};

	double GetItemWeight(Item *I) 
	{
		return I->Weight;
	};


	struct Cluster 
	{
		// Value of the energy term for this cluster (cached to increase speed by about 33%)
		double EnergyValue;

		// The accumulated Normals (Weighted)
		double SNormal[3];

		// The accumulated centroids (weighted)
		double SValue[3];

		// the accumulated weights
		double SWeight;

	};

	double GetClusterEnergy(Cluster *C)
	{
		C->EnergyValue=
            (-C->SNormal[0]*C->SNormal[0]
			-C->SNormal[1]*C->SNormal[1]
			-C->SNormal[2]*C->SNormal[2]
			-Factor*
				(
				C->SValue[0]*C->SValue[0]
				+C->SValue[1]*C->SValue[1]
				+C->SValue[2]*C->SValue[2]
				)
				)/C->SWeight;

		return C->EnergyValue;
	}


	int GetClusterRankDeficiency(Cluster *C)
	{
		return (0);
	}
	void ComputeClusterEnergy(Cluster *C)
	{
	}

	double ComputeDistanceBetweenItemAndCluster(Item *I,Cluster *C)
	{
		double P1[3],P2[3];
		int i;
		for (i=0;i<3;i++)
		{
			P1[i]=I->Normal[i];
			P2[i]=C->SNormal[i];
		}

		vtkMath::Normalize(P1);
		vtkMath::Normalize(P2);
		return vtkMath::Distance2BetweenPoints(P1,P2);

	}

	double ComputeDistanceBetweenItemAndPoint(Item *I,double *P1)
	{
		cout<<"should never happen! L21MetricforClustering"<<endl;

		double P2[3];
		P2[0]=I->Normal[0]/I->Weight;
		P2[1]=I->Normal[1]/I->Weight;
		P2[2]=I->Normal[2]/I->Weight;
		vtkMath::Normalize(P2);

		return vtkMath::Distance2BetweenPoints(P1,P2);
	}

	void DeepCopy(Cluster *Source,Cluster *Destination)
	{
		int i;
		Destination->SWeight=Source->SWeight;
		for (i=0;i<3;i++)
		{
			Destination->SValue[i]=Source->SValue[i];
			Destination->SNormal[i]=Source->SNormal[i];
		}

		Destination->EnergyValue=Source->EnergyValue;
	}

	void AddItemToCluster(Item *I,Cluster *C)
	{
		int i;
		for (i=0;i<3;i++)
		{
			C->SValue[i]+=I->Value[i];
			C->SNormal[i]+=I->Normal[i];
		}
		C->SWeight+=I->Weight;
	}

	void SubstractItemFromCluster(Item *I,Cluster *C)
	{
		int i;
		for (i=0;i<3;i++)
		{
			C->SValue[i]-=I->Value[i];
			C->SNormal[i]-=I->Normal[i];
		}
		C->SWeight-=I->Weight;
	}

	void ComputeClusterCentroid(Cluster *C)
	{
	}

	void GetClusterCentroid(Cluster *C,double	*P)
	{
		P[0]=C->SValue[0]/C->SWeight;
		P[1]=C->SValue[1]/C->SWeight;
		P[2]=C->SValue[2]/C->SWeight;
	}

	void ResetCluster(Cluster *C) 
	{
		int i;
		for (i=0;i<3;i++)
		{
			C->SValue[i]=0;
			C->SNormal[i]=0;
		}
		C->SWeight=0;
		C->EnergyValue=0;
	}
	vtkL21MetricForClustering()
	{
		this->CustomWeights=0;
		this->Gradation=0;
	};
	~vtkL21MetricForClustering()
	{
	};

	void BuildMetric( vtkSurface *Mesh, int ClusteringType )
	{
		vtkIdType i,j;
		vtkIdType v1,v2,v3;
		double P1[3],P2[3],P3[3];

		// Compute the parameter mixing between the L21 metric and the isotropic metric
		// We actually mix L21 with Isotropic metric to improve convergence
		// The mixing parameter is chosen so that the Isotropic metric is negligible
		// except in flat regions

		double Bounds[6];
		Mesh->ComputeBounds();
		Mesh->GetBounds(Bounds);
		P1[0]=Bounds[0];
		P2[0]=Bounds[1];
		P1[1]=Bounds[2];
		P2[1]=Bounds[3];
		P1[2]=Bounds[4];
		P2[2]=Bounds[5];

		Factor=0;//0.0001/vtkMath::Distance2BetweenPoints(P1,P2);


		if (ClusteringType==0)
		{
			// Items are triangles
			Items=new Item[Mesh->GetNumberOfCells()];
			for (i=0;i<Mesh->GetNumberOfCells();i++)
			{
				this->ResetItem(Items+i);
				Mesh->GetFaceVertices(i,v1,v2,v3);
				Mesh->GetPoint(v1,P1);
				Mesh->GetPoint(v2,P2);
				Mesh->GetPoint(v3,P3);
				if (CustomWeights==0)
					Items[i].Weight=vtkTriangle::TriangleArea(P1,P2,P3);
				else
					Items[i].Weight=vtkTriangle::TriangleArea(P1,P2,P3)*pow(CustomWeights->GetValue(i),Gradation);
				Items[i].Value[0]=(P1[0]+P2[0]+P3[0])/3.0;
				Items[i].Value[1]=(P1[1]+P2[1]+P3[1])/3.0;
				Items[i].Value[2]=(P1[2]+P2[2]+P3[2])/3.0;

				vtkTriangle::ComputeNormal(P1,P2,P3,Items[i].Normal);
			}
			this->ClampWeights(Items,Mesh->GetNumberOfCells(),10000);
			for (i=0;i<Mesh->GetNumberOfCells();i++)
			{
				for (j=0;j<3;j++)
				{
					Items[i].Value[j]*=Items[i].Weight;
					Items[i].Normal[j]*=Items[i].Weight;
				}
			}
		}
		else
		{
			// Items are vertices
			vtkIdList *FList=vtkIdList::New();
			double N[3];

			Items=new Item[Mesh->GetNumberOfPoints()];
			vtkDoubleArray *VerticesAreas=Mesh->GetVerticesAreas();
			for (i=0;i<Mesh->GetNumberOfPoints();i++)
			{
				this->ResetItem(Items+i);

				if (CustomWeights==0)
					Items[i].Weight=VerticesAreas->GetValue(i);
				else
					Items[i].Weight=VerticesAreas->GetValue(i)*pow(CustomWeights->GetValue(i),Gradation);

				Mesh->GetVertexNeighbourFaces(i,FList);
				Mesh->GetPointCoordinates(i,Items[i].Value);
				for (j=0;j<FList->GetNumberOfIds();j++)
				{
					Mesh->GetFaceVertices(FList->GetId(j),v1,v2,v3);
					Mesh->GetPoint(v1,P1);
					Mesh->GetPoint(v2,P2);
					Mesh->GetPoint(v3,P3);
					double TriangleArea=vtkTriangle::TriangleArea(P1,P2,P3);
					vtkTriangle::ComputeNormalDirection(P1,P2,P3,N);
					Items[i].Normal[0]+=N[0]*TriangleArea;
					Items[i].Normal[1]+=N[1]*TriangleArea;
					Items[i].Normal[2]+=N[2]*TriangleArea;
				}
				vtkMath::Normalize(Items[i].Normal);
			}
			FList->Delete();
			this->ClampWeights(Items,Mesh->GetNumberOfPoints(),10000);
			for (i=0;i<Mesh->GetNumberOfPoints();i++)
			{
				for (j=0;j<3;j++)
				{
					Items[i].Value[j]*=Items[i].Weight;
					Items[i].Normal[j]*=Items[i].Weight;
				}
			}
		}
	}
	// this method clamps the weights between AverageValue/Ratio and AverageValue*Ratio
	void ClampWeights(Item *Items,int NumberOfValues,double Ratio)
	{
		double Average=0;
		int i;
		for (i=0;i<NumberOfValues;i++)
		{
			Average+=Items[i].Weight;
		}
		Average=Average/(double)NumberOfValues;
		double Min=Average/Ratio;
		double Max=Average*Ratio;
		for (i=0;i<NumberOfValues;i++)
		{
			if (Items[i].Weight>Max)
				Items[i].Weight=Max;
			if (Items[i].Weight<Min)
				Items[i].Weight=Min;
		}
	}
private:
	vtkDoubleArray *CustomWeights;
	double Factor;
	double Gradation;

	// The array storing items
	Item *Items;

};

#endif
