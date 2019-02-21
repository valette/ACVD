/*=========================================================================

Program:   Lloyd relaxation based Clustering for meshes 
Module:    vtkLloydClustering.h
Language:  C++
Date:      2007/03
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

#ifndef _VTKTLLOYDCLUSTERING_H_
#define _VTKTLLOYDCLUSTERING_H_

#include <vtkPriorityQueue.h>
#include "vtkUniformClustering.h"

class STDPriorityQueue;
//typedef vtkPriorityQueue MyPriorityQueue;
typedef STDPriorityQueue MyPriorityQueue;

template < class Metric > class vtkLloydClustering:public vtkUniformClustering <	Metric >
{
protected:

	virtual int ProcessOneLoop();
	virtual void Init();

	void PushItemRing(int Item);

	vtkIntArray *ItemsLastVisitTime;
	vtkIntArray *ClustersClosestItem;
	vtkDoubleArray *ClustersClosestItemDistance;


	class LloydPriorityQueue
	{
	public:

		void Allocate(int N) {};
		int GetNumberOfItems(){return Queue.size();};
		void Insert(double Priority,int Item)
		{
			Distance Dist;
			Dist.Item=Item;
			Dist.Priority=Priority;
			Queue.push(Dist);
		}
		
		int Pop()
		{
			Distance Dist;
			Dist=Queue.top();
			Queue.pop();
			return (Dist.Item);
		}
		class Distance
		{
		public:
			double Priority;
			vtkIdType Item;
			
			bool operator < (const Distance& d) const
		   {
		      return (Priority > d.Priority);
		   }			
		};

		static LloydPriorityQueue* New() {return new LloydPriorityQueue;};
		void Delete() {delete (this);};
		std::priority_queue < Distance> Queue;

		// the constructor
		LloydPriorityQueue (){};

		// the destructor
		~LloydPriorityQueue (){};
	};

	LloydPriorityQueue *PriorityQueue;

	// the constructor
	vtkLloydClustering ();
	
	// the destructor
	~vtkLloydClustering ();
};

template < class Metric >
void	vtkLloydClustering < Metric >::PushItemRing (int Item)
{
	int I1,I2;
	int Edge;
	double P1[3],P2[3];
	
	this->MetricContext.GetClusterCentroid(this->Clusters+this->Clustering->GetValue(Item),P1);
	vtkIdList *EList=vtkIdList::New();
	this->GetItemEdges(Item,EList);
	
	int i;
	for (i=0;i<EList->GetNumberOfIds();i++)
	{
		Edge=EList->GetId(i);
		this->GetEdgeItems(Edge,I1,I2);
		if (I2==Item)
		{
			I2=I1;
		}
		if (this->ItemsLastVisitTime->GetValue(I2)<this->NumberOfLoops)
		{
			// the Item I2 was not previously visited so we push the edge in the queue
			this->GetItemCoordinates(I2,P2);
			this->PriorityQueue->Insert(vtkMath::Distance2BetweenPoints(P1,P2),Edge);
		}
	}	
	EList->Delete();
}

template < class Metric >
int	vtkLloydClustering < Metric >::ProcessOneLoop ()
{	
	int Item;
	int Cluster;
	double Distance2;
	int Edge;
	int NumberOfModifications=0;
	
	double P1[3],P2[3];
	
	// first compute for each cluster the closest Item;
	for (Cluster=0;Cluster<this->NumberOfClusters;Cluster++)
	{
		this->ClustersClosestItemDistance->SetValue(Cluster,-1);
	}

	for (Item=0;Item<this->GetNumberOfItems();Item++)
	{
		this->ItemsLastVisitTime->SetValue(Item,this->NumberOfLoops-1);
		Cluster=this->Clustering->GetValue(Item);
		if (Cluster<this->NumberOfClusters)
		{
			this->MetricContext.GetClusterCentroid(this->Clusters+Cluster,P1);
			this->GetItemCoordinates(Item,P2);
			Distance2=vtkMath::Distance2BetweenPoints(P1,P2);
			if ((this->ClustersClosestItemDistance->GetValue(Cluster)<0)
				||(Distance2<this->ClustersClosestItemDistance->GetValue(Cluster)))
			{
				this->ClustersClosestItemDistance->SetValue(Cluster,Distance2);
				this->ClustersClosestItem->SetValue(Cluster,Item);
			}
		}
	}
	
	vtkTimerLog *Timer=vtkTimerLog::New();
	Timer->StartTimer();
	
	// push for each cluster the ring of each closest item
	for (Cluster=0;Cluster<this->NumberOfClusters;Cluster++)
	{
		Item=this->ClustersClosestItem->GetValue(Cluster);
		this->ItemsLastVisitTime->SetValue(Item,this->NumberOfLoops);
		this->PushItemRing(Item);
	}
		
	// create the Voronoi diagram with the priority queue	
	int TempItem;
	while (this->PriorityQueue->GetNumberOfItems())
	{
		int I1,I2;
		Edge=this->PriorityQueue->Pop();
		this->GetEdgeItems(Edge,I1,I2);
		if (this->ItemsLastVisitTime->GetValue(I2)==this->NumberOfLoops)
		{
			TempItem=I1;
			I1=I2;
			I2=TempItem;
		}
		if (this->ItemsLastVisitTime->GetValue(I2)<this->NumberOfLoops)
		{
			Cluster=this->Clustering->GetValue(I1);
			if (this->Clustering->GetValue(I2)!=Cluster)
				NumberOfModifications++;
			this->Clustering->SetValue(I2,Cluster);
			this->ItemsLastVisitTime->SetValue(I2,this->NumberOfLoops);
			this->PushItemRing(I2);
		}
	}
	
	Timer->StopTimer();
	if (this->ConsoleOutput)
	{
		cout<<endl<<Timer->GetElapsedTime()<<" seconds for priority queue handling "<<endl;
	}
	this->ReComputeStatistics();
	return (NumberOfModifications);
}

template < class Metric > void vtkLloydClustering < Metric >::Init ()
{
	vtkUniformClustering < Metric >::Init ();
	this->ClustersClosestItemDistance->SetNumberOfValues(this->NumberOfClusters);
	this->ClustersClosestItem->SetNumberOfValues(this->NumberOfClusters);	
	this->ItemsLastVisitTime->SetNumberOfValues(this->GetNumberOfItems());
	
	int Item;
	for (Item=0;Item<this->GetNumberOfItems();Item++)
	{
		this->ItemsLastVisitTime->SetValue(Item,this->NumberOfLoops-1);
	}
	this->PriorityQueue->Allocate(2*this->GetNumberOfEdges());
}


template < class Metric >
	vtkLloydClustering < Metric >::vtkLloydClustering ()
{
	this->ItemsLastVisitTime=vtkIntArray::New();
	this->ClustersClosestItem=vtkIntArray::New();
	this->ClustersClosestItemDistance=vtkDoubleArray::New();
	this->PriorityQueue=LloydPriorityQueue::New();
}

template < class Metric >
	vtkLloydClustering < Metric >::~vtkLloydClustering ()
{
	this->ItemsLastVisitTime->Delete();
	this->ClustersClosestItem->Delete();
	this->ClustersClosestItemDistance->Delete();
	this->PriorityQueue->Delete();
}

#endif
