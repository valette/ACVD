/***************************************************************************
                          vtkClusteringExampleFile.h  -  description
                             -------------------
    begin                : October 2006
    copyright            : (C) 2006 by Sebastien Valette
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
#ifndef _VTKCLUSTERINGEXAMPLEFILE_H_
#define _VTKCLUSTERINGEXAMPLEFILE_H_

#include <vtkObjectFactory.h>
#include "vtkUniformClustering.h"


template < class Metric, typename EdgeType=vtkIdType > class
vtkClusteringExample:public vtkUniformClustering <Metric,EdgeType>
{

public:
	static vtkClusteringExample *New()
	{
		vtkObject* ret = vtkObjectFactory::CreateInstance("vtkClusteringExample");
		if(ret)
	    {
			return (vtkClusteringExample*)ret;
	    }
	  // If the factory was unable to create the object, then create it here.
		return (new vtkClusteringExample<Metric>);
	};


protected:
	vtkClusteringExample() 
	{
	};

	~vtkClusteringExample()
	{};

	void GetItemCoordinates(vtkIdType Item, double *P)
	{
	}
	
	double GetItemArea(vtkIdType Item)
	{
	return (0);
	}

	int GetNumberOfDualItems()
	{return (0);};




	void GetItemNeighbours (vtkIdType Item, vtkIdList *IList)
	{
	};

	int GetNumberOfItems ()
	{
	return (0);
	};

	void AddItemRingToProcess(vtkIdType Item)
	{
//			this->AddEdgeToProcess();
	};

	void GetEdgeItems(vtkIdType Item, vtkIdType &I1,vtkIdType &I2)
	{
	};

	void GetItemEdges(vtkIdType Item, vtkIdList *EList)
	{
	}

	/// check if moving the face Item from Cluster will not make Cluster non-connex
	/// returns 0 if no problem; 1 otherwise
	int ConnexityConstraintProblem(vtkIdType Item,EdgeType Edge,vtkIdType Cluster, vtkIdType Cluster2)
	{
		if (!this->ConnexityConstraint)
			return (0);
	
	return (0);
	};


	EdgeType GetNumberOfEdges()
	{
	return (0);
	};

	void GetEdgeItemsSure (vtkIdType Item, vtkIdList * VList)
	{};
	

	void SetNumberOfProcesses (int N)
	{};

protected:

	int NumberOfProcesses;

	void CreateWindows(){};
	void Snapshot(){};
	void BuildMetric(){};
};

class  vtkMetricExample
{
public:

	void SetConstrainedClustering(int C){};

	struct Item
	{
		double Weight;
	};

	void ResetItem(Item *I) 
	{
	}

	double GetItemWeight(vtkIdType Item) 
	{
		return 0;
	}

	struct Cluster 
	{
		double SWeight;
	};

	double GetClusterEnergy(Cluster *C)
	{
	return (0);
	}

	void ComputeClusterEnergy(Cluster *C)
	{
	}

	void DeepCopy(Cluster *Source, Cluster *Destination)
	{
	}

	void Add(Cluster *Source, vtkIdType ItemId, Cluster *Destination)
	{
	}

	void AddItemToCluster(vtkIdType Item, Cluster *C)
	{
	}

	void Sub(Cluster *Source, vtkIdType ItemId, Cluster *Destination)
	{
	}

	void SubstractItemFromCluster(vtkIdType Item,Cluster *C)
	{
	}
	void ComputeClusterCentroid(Cluster *C)
	{
	}

	void GetClusterCentroid(Cluster *C, double *P)
	{
	}
	void SetClusterCentroid(Cluster *C, double *P)
	{
	}
	void ResetCluster(Cluster *C) 
	{
	}

	vtkMetricExample()
	{
	};
	~vtkMetricExample()
	{
	};
};

#endif

