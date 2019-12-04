/*=========================================================================

Program:   Isotropic Metric for Clustering 
Module:    vtkIsotropicMetricForClustering.h
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

#ifndef _VTKISOTROPICMETRICFORCLUSTERING_H_
#define _VTKISOTROPICMETRICFORCLUSTERING_H_

#include <vtkTriangle.h>
#include <vtkDataArrayCollection.h>
#include "vtkSurface.h"

class  vtkIsotropicMetricForClustering
{
public:

	// this method is empty as there is no constraint for this metric
	void SetConstrainedClustering(int C){}

	void MultiplyItemWeight(vtkIdType ItemId, double Factor)
	{
		for (int i=0;i<3;i++)
			this->Items[ItemId].Value[i]*=Factor;
		this->Items[ItemId].Weight*=Factor;
	}

	int IsCurvatureIndicatorNeeded()
	{
		if (this->Gradation>0)
			return 1;
		else
			return (0);
	}
	int IsPrincipalDirectionsNeeded()
	{
		return (0);
	}

	void SetCurvatureInfo(vtkDataArrayCollection *Info)
	{
		if (this->CustomWeights)
			this->CustomWeights->Delete();
		this->CustomWeights=(vtkDoubleArray*) Info->GetItem(0);
		this->CustomWeights->Register(this->Object);
	}
	void SetGradation (double Gradation)
	{
		this->Gradation=Gradation;
	}

	double GetGradation ()
	{
		return this->Gradation;
	}

	struct Item
	{
		double Value[3];
		double Weight;
	};

	void ResetItem(Item *I) 
	{
		I->Value[0]=0;
		I->Value[1]=0;
		I->Value[2]=0;
		I->Weight=0;
	};
	double GetItemWeight(vtkIdType ItemId) 
	{
		return this->Items[ItemId].Weight;
	};

	struct Cluster 
	{
		double SValue[3];
		double SWeight;
		double EnergyValue;
	};

	int GetClusterRankDeficiency(Cluster *C)
	{return (0);};
	double GetClusterEnergy(Cluster *C)
	{
		return C->EnergyValue;
	}
	void ComputeClusterEnergy(Cluster *C)
	{
		C->EnergyValue=(
			-C->SValue[0]*C->SValue[0]
			-C->SValue[1]*C->SValue[1]
			-C->SValue[2]*C->SValue[2]
			)/C->SWeight;
	}
	double ComputeDistanceBetweenItemAndCluster(Item *I,Cluster *C)
	{
		return (0);
	}

	double ComputeDistanceBetweenItemAndPoint(Item *I,double *P1)
	{
		double P2[3];
		P2[0]=I->Value[0]/I->Weight;
		P2[1]=I->Value[1]/I->Weight;
		P2[2]=I->Value[2]/I->Weight;
		return sqrt(vtkMath::Distance2BetweenPoints(P1,P2));
	}


	void DeepCopy(Cluster *Source, Cluster *Destination)
	{
		Destination->SWeight=Source->SWeight;
		Destination->SValue[0]=Source->SValue[0];
		Destination->SValue[1]=Source->SValue[1];
		Destination->SValue[2]=Source->SValue[2];
		Destination->EnergyValue=Source->EnergyValue;
	}

	void Add(Cluster *Source, vtkIdType ItemId, Cluster *Destination)
	{
		this->DeepCopy(Source,Destination);
		this->AddItemToCluster(ItemId, Destination);
	}

	void AddItemToCluster(vtkIdType ItemId, Cluster *C)
	{
		Item *I=&this->Items[ItemId];
		C->SValue[0]+=I->Value[0];
		C->SValue[1]+=I->Value[1];
		C->SValue[2]+=I->Value[2];
		C->SWeight+=I->Weight;
	}

	void Sub(Cluster *Source, vtkIdType ItemId, Cluster *Destination)
	{
		this->DeepCopy(Source,Destination);
		this->SubstractItemFromCluster(ItemId, Destination);		
	}

	void SubstractItemFromCluster(vtkIdType ItemId, Cluster *C)
	{
		Item *I=&this->Items[ItemId];
		C->SValue[0]-=I->Value[0];
		C->SValue[1]-=I->Value[1];
		C->SValue[2]-=I->Value[2];
		C->SWeight-=I->Weight;
	}
	void ComputeClusterCentroid(Cluster *C)
	{}

	void GetClusterCentroid(Cluster *C, double *P)
	{
		P[0]=C->SValue[0]/C->SWeight;
		P[1]=C->SValue[1]/C->SWeight;
		P[2]=C->SValue[2]/C->SWeight;
	}
	void SetClusterCentroid(Cluster *C, double *P)
	{
		C->SValue[0]=P[0]*C->SWeight;
		C->SValue[1]=P[1]*C->SWeight;
		C->SValue[2]=P[2]*C->SWeight;
	}
	void ResetCluster(Cluster *C) 
	{
		C->SValue[0]=0; 
		C->SValue[1]=0; 
		C->SValue[2]=0;
		C->SWeight=0;
		C->EnergyValue=0;
	}
	vtkIsotropicMetricForClustering()
	{
		this->CustomWeights=0;
		this->Gradation=0;
		this->Object=vtkObject::New();
		this->Items=0;
	};

	~vtkIsotropicMetricForClustering()
	{
		this->Object->Delete();
		if (this->CustomWeights)
			this->CustomWeights->Delete();
		if (this->Items)
			delete [] this->Items;
	};

	void BuildMetric( vtkSurface *Mesh, int ClusteringType )
	{
		vtkIdType i;

		//Build the items
		if (ClusteringType==0)
		{
			// Items are triangles
			vtkIdType v1,v2,v3;
			double P1[3],P2[3],P3[3];
			Items=new Item[Mesh->GetNumberOfCells()];
			for (i=0;i<Mesh->GetNumberOfCells();i++)
			{
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

			}
			this->ClampWeights(Items,Mesh->GetNumberOfCells(),100000);
			for (i=0;i<Mesh->GetNumberOfCells();i++)
			{
				Items[i].Value[0]*=Items[i].Weight;
				Items[i].Value[1]*=Items[i].Weight;
				Items[i].Value[2]*=Items[i].Weight;
			}
		}
		else
		{
			// Items are vertices
			Items=new Item[Mesh->GetNumberOfPoints()];
			vtkDoubleArray *VerticesAreas=Mesh->GetVerticesAreas();
			for (i=0;i<Mesh->GetNumberOfPoints();i++)
			{
				if (CustomWeights==0)
					Items[i].Weight=VerticesAreas->GetValue(i);
				else
					Items[i].Weight=VerticesAreas->GetValue(i)*pow(CustomWeights->GetValue(i),Gradation);

				Mesh->GetPoint(i,Items[i].Value);
			}
			this->ClampWeights(Items,Mesh->GetNumberOfPoints(),100000);
			for (i=0;i<Mesh->GetNumberOfPoints();i++)
			{
				Items[i].Value[0]*=Items[i].Weight;
				Items[i].Value[1]*=Items[i].Weight;
				Items[i].Value[2]*=Items[i].Weight;
			}
		}
	}
	// this method clamps the weights between AverageValue/Ratio and AverageValue*Ratio
	static void ClampWeights(Item *Items,int NumberOfValues,double Ratio)
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
	
	// Dummy object used for registering the curvature indicators
	vtkObject	*Object;
	double Gradation;
	
	// The array storing items
	Item *Items;
};

#endif
