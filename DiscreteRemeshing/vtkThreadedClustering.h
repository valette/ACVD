/*=========================================================================

Program:   Uniform Clustering for meshes (abstract class, multithreaded version)
Module:    vtkThreadedClustering.h
Language:  C++
Date:      2006/01
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

// The template parameter Derived (the most-derived
// clustering class, CRTP) makes the per-edge calls of ExecuteProcess and
// ExecuteProcessWithDistances static instead of virtual.
// Each loop runs the regions in conflict-free rounds
// (vtkSMPTools) instead of a vtkMultiThreader pool ordered by measured run
// times with per-cluster mutexes, so the clustering is repeatable and does not
// depend on the number of threads. The region count follows from
// ItemsPerRegion and the cluster size, not from the number of threads.

#ifndef _VTKTHREADEDCLUSTERING_H_
#define _VTKTHREADEDCLUSTERING_H_

#include <vtkMath.h>
#include <vtkSMPTools.h>
#include <vtkTimerLog.h>
#include <algorithm>
#include <deque>
#include <vector>

#include "vtkUniformClustering.h"

// Class derived from vtkUniformClustering

template < class Metric, class Derived > class vtkThreadedClustering:public vtkUniformClustering <	Metric >
{
public:

	/// Sets the least number of items of a region (one parallel job). A region
	/// is also at least RegionClusterWidth clusters large. The result depends
	/// on it; the number of threads only changes the speed.
	void SetItemsPerRegion (int N)
	{
		this->ItemsPerRegion=N;
	};

	int GetItemsPerRegion ()
	{
		return (this->ItemsPerRegion);
	};

protected:

	// the most-derived object, for the statically dispatched per-edge calls
	Derived* Self() { return static_cast<Derived*>(this); }

	// the constructor
	vtkThreadedClustering ();

	// the destructor
	~vtkThreadedClustering ();

	// Executes one process (a region, or the edges between regions)
	void ExecuteProcess (int Process);
	void ExecuteProcessWithDistances (int Process);
	void ExecuteRegion (int Process);

	// collects in RegionClusters[Region] the clusters Region can touch in this loop
	void CollectRegionClusters (int Region);

	// splits the regions into Rounds whose regions touch disjoint clusters and are not adjacent
	void ComputeRounds ();

	// virtual function. Might be implemented in derived classes for speed issues.
	// Thread is the scratch slot, one per process.
	virtual void AddItemRingToProcess (vtkIdType Item, int ProcessId, int Thread)
	{
		vtkIdList *EList=this->ProcessInfos[Thread].itemList;
		this->Self()->Derived::GetItemEdges(Item,EList);
		for (int i=0;i<EList->GetNumberOfIds();i++)
			this->vtkThreadedClustering::AddEdgeToProcess (EList->GetId(i),ProcessId);
	};

	// virtual function. Might be implemented in derived classes for speed issues
	virtual int ThreadedConnexityConstraintProblem(vtkIdType Item,vtkIdType Edge,vtkIdType Cluster2, vtkIdType Cluster1 ,int Thread)
		{return this->ConnexityConstraintProblem(Item,Edge,Cluster2,Cluster1);};

	// allocates and initializes memory for the threaded clustering
	// (queues, timings)
	virtual void Init ();

	// returns the Items adjacent to the given edge ("Sure" means that
	// the non-manifold cases are well managed)
	virtual void GetEdgeItemsSure (vtkIdType Item, vtkIdList * VList)=0;

	// assigns the edges to processes in the threads pool. The result is stored in the
	// EdgesProcess Array
	void ComputeEdgesLayout ();

	// a sub-part of edges layout computation
	void ComputeMeshSlices (int NumberOfSlices,vtkIntArray * Clustering);

	/// this method replaces the AddEdgeToProcess in a multithreaded context
	virtual void AddEdgeToProcess (vtkIdType Edge, int ProcessId)
	{
		this->ProcessesPushQueues[this->
					  EdgesProcess[Edge]][ProcessId].push_back (Edge);
	};

	/// The number of processes: the regions plus the one for the edges between them
	int PoolSize;

	/// The least number of items of a region
	int ItemsPerRegion;

	/// The least size of a region in clusters, so that clusters rarely span non-adjacent regions
	static constexpr int RegionClusterWidth = 8;

	/// The method used to compute the edges layout
	int EdgesLayoutComputingType;

	/// The main funcion : processes one loop on the boundary edges
	virtual int ProcessOneLoop();

	// parameter to set On/Off the threads timings (0:Off 1: On)
	int DisplayThreadsTimingsFlag;

	/// Method which displays on screen the execution times for each thread
	void DisplayThreadsTimings();

	// refills the queues according to a possibly updated clustering
	void FillQueuesFromClustering ();

	// Context for Clustering
	// *******************************************
	std::vector<int> EdgesProcess;

	// the regions sharing an edge with each region (sorted)
	std::vector< std::vector<int> > RegionNeighbours;

	class ProcessData {
		public:
		int PreviousNumberOfIterations;
		int NumberOfIterations;
		int NumberOfModifications;
		double StartTime;
		double StopTime;
		vtkIdList *itemList;
		ProcessData() { itemList = vtkIdList::New(); };
		~ProcessData() { itemList->Delete(); };
	};

	std::vector< ProcessData > ProcessInfos;

	// Schedule of one loop
	// **********************************************
	std::vector< std::vector<vtkIdType> > RegionClusters;
	std::vector<int> ClusterRound;
	std::vector<int> RegionRound;
	std::vector< std::vector<int> > Rounds;

	// The push and pop queues used to track the edges laying between different clusters.
	// Queue [target process][pushing process]: each process pushes into its own column.
	std::deque <int>**ProcessesQueues1;
	std::deque <int>**ProcessesQueues2;

	std::deque <int>**ProcessesPushQueues;
	std::deque <int>**ProcessesPopQueues;

	// this method swaps the pop queues with the push queues.
	void SwapQueues ();
};

template < class Metric, class Derived > void vtkThreadedClustering < Metric, Derived >::FillQueuesFromClustering ()
{
	int i, j;
	vtkIdType I1, I2;

	// empty the pushing queues
	for (i = 0; i < this->PoolSize; i++)
	{
		for (j = 0; j < this->PoolSize; j++)
		{
			this->ProcessesPushQueues[i][j].clear ();
		}
	}

	// fill the pushqueues depending on the clustering
	for (i = 0; i < this->GetNumberOfEdges (); i++)
	{
		this->GetEdgeItems (i, I1, I2);
		if (I2 >= 0)
		{
			if (this->Clustering->GetValue (I1) !=this->Clustering->GetValue (I2))
			{
				this->ProcessesPushQueues[this->EdgesProcess[i]]
					[this->EdgesProcess[i]].push_back (i);
			}
		}
	}
}

// A region writes only its own items, and moves an item only into the cluster
// of the other end of an edge it processes. So the clusters it touches in a
// loop are those of the ends of its queued edges whose clusters differ at loop
// start (an edge that is uniform at start can differ later only after one of
// its ends moved, which touched that cluster already).
template < class Metric, class Derived > void vtkThreadedClustering < Metric, Derived >::CollectRegionClusters (int Region)
{
	std::vector<vtkIdType> &Clusters = this->RegionClusters[Region];
	Clusters.clear ();
	for (int CurrentQueue = 0; CurrentQueue < this->PoolSize; CurrentQueue++)
	{
		for (int Edge : this->ProcessesPopQueues[Region][CurrentQueue])
		{
			vtkIdType I1, I2;
			this->Self()->Derived::GetEdgeItems (Edge, I1, I2);
			if (I2 < 0) continue;
			vtkIdType Val1 = this->Clustering->GetValue (I1);
			vtkIdType Val2 = this->Clustering->GetValue (I2);
			if (Val1 == Val2) continue;
			// NumberOfClusters marks an unassigned item, it is no cluster state
			if (Val1 != this->NumberOfClusters) Clusters.push_back (Val1);
			if (Val2 != this->NumberOfClusters) Clusters.push_back (Val2);
		}
	}
}

// Greedy in region order: a region joins the first round in which it shares no
// cluster and no edge with a region already there. Adjacent regions are kept
// apart because the connexity test reads the clusters of an item's ring.
template < class Metric, class Derived > void vtkThreadedClustering < Metric, Derived >::ComputeRounds ()
{
	const int NumberOfRegions = this->PoolSize - 1;
	this->ClusterRound.assign (this->NumberOfClusters, -1);
	this->RegionRound.assign (NumberOfRegions, -1);
	this->Rounds.clear ();

	std::vector<int> Pending (NumberOfRegions), Deferred;
	for (int Region = 0; Region < NumberOfRegions; Region++) Pending[Region] = Region;

	while (!Pending.empty ())
	{
		const int Round = static_cast<int> (this->Rounds.size ());
		this->Rounds.emplace_back ();
		Deferred.clear ();
		for (int Region : Pending)
		{
			bool Free = true;
			for (int Neighbour : this->RegionNeighbours[Region])
				if (this->RegionRound[Neighbour] == Round) { Free = false; break; }
			if (Free)
				for (vtkIdType Cluster : this->RegionClusters[Region])
					if (this->ClusterRound[Cluster] == Round) { Free = false; break; }
			if (!Free)
			{
				Deferred.push_back (Region);
				continue;
			}
			this->RegionRound[Region] = Round;
			for (vtkIdType Cluster : this->RegionClusters[Region])
				this->ClusterRound[Cluster] = Round;
			this->Rounds.back ().push_back (Region);
		}
		Pending.swap (Deferred);
	}
}

template < class Metric, class Derived > void vtkThreadedClustering < Metric, Derived >::ExecuteRegion (int Process)
{
	ProcessData &data = this->ProcessInfos[ Process ];
	data.StartTime = vtkTimerLog::GetUniversalTime ();
	if (this->MinimizeUsingEnergy)
		this->ExecuteProcess (Process);
	else
		this->ExecuteProcessWithDistances (Process);
	data.StopTime = vtkTimerLog::GetUniversalTime ();
}

template < class Metric, class Derived > void
	vtkThreadedClustering < Metric, Derived >::ExecuteProcess (int Process)
{
	vtkIdType Edge, I1, I2;
	int Val1, Val2, *Size1, *Size2;

	// Those variables will contain the energy values for the three possible cases.
	// The "volatile" statement is here to fix some numerical issues
	// (see : http://gcc.gnu.org/bugzilla/show_bug.cgi?id=323)
	volatile double Try1, Try2, Try3;
	volatile double Try11, Try21, Try31;
	volatile double Try12, Try22, Try32;

	typename Metric::Cluster *Cluster1, *Cluster2, *Cluster21, *Cluster22, *Cluster31,*Cluster32;
	Cluster21 = new typename Metric::Cluster;
	Cluster22 = new typename Metric::Cluster;
	Cluster31 = new typename Metric::Cluster;
	Cluster32 = new typename Metric::Cluster;
	ProcessData &data = this->ProcessInfos[ Process ];

	std::deque < int > *Queue;
	for (int CurrentQueue=0;CurrentQueue<this->PoolSize;CurrentQueue++)
	{
		Queue = &this->ProcessesPopQueues[Process][CurrentQueue];
		while (!Queue->empty())
		{
			// the queue is not empty : let's use it
			Edge = Queue->front ();
			Queue->pop_front ();
			this->Self()->Derived::GetEdgeItems (Edge, I1, I2);

			if ((this->EdgesLastLoop[Edge] ==this->RelativeNumberOfLoops) || (I2 < 0)) continue;
			this->EdgesLastLoop[Edge] = this->RelativeNumberOfLoops;
			Val1 = this->Clustering-> GetValue (I1);
			Val2 = this->Clustering-> GetValue (I2);
			if (Val2==Val1) continue;

			data.NumberOfIterations++;
			if (Val1 == this->NumberOfClusters)
			{
				// I1 is not associated. Give it to the same cluster as I2
				this->MetricContext.AddItemToCluster(I1,&this->Clusters[Val2]);
				this->MetricContext.ComputeClusterCentroid(&this->Clusters[Val2]);
				this->MetricContext.ComputeClusterEnergy(&this->Clusters[Val2]);
				(*this->ClustersSizes->GetPointer (Val2))++;
				this->vtkThreadedClustering::AddItemRingToProcess (I1,Process, Process);
				this->Clustering->SetValue (I1,Val2);
				data.NumberOfModifications++;
			}
			else if (Val2 ==this->NumberOfClusters)
			{
				// I2 is not associated. Give it to the same cluster as I1
				this->MetricContext.AddItemToCluster(I2,&this->Clusters[Val1]);
				this->MetricContext.ComputeClusterCentroid(&this->Clusters[Val1]);
				this->MetricContext.ComputeClusterEnergy(&this->Clusters[Val1]);
				(*this->ClustersSizes->GetPointer (Val1))++;
				this->vtkThreadedClustering::AddItemRingToProcess (I2,Process,Process);
				this->Clustering->SetValue (I2,Val1);
				data.NumberOfModifications++;
			}
			else
			{
				int Result;

				// determine whether one of	the	two	adjacent clusters was modified,
				// or whether any of the clusters is freezed
				//	If not,	the	test is	useless, and the speed improved	:)
				if (((this->ClustersLastModification[Val1] >=this->NumberOfLoops-1)
					|| (this->ClustersLastModification[Val2]>=this->NumberOfLoops-1))
					&&((this->IsClusterFreezed->GetValue(Val1)==0)
						&&(this->IsClusterFreezed->GetValue(Val2)==0)))
				{
					Cluster1 = &this->Clusters[Val1];
					Cluster2 = &this->Clusters[Val2];

					Size1 = this->ClustersSizes->GetPointer (Val1);
					Size2 = this->ClustersSizes->GetPointer (Val2);

					// Compute the initial energy
					Try11=this->MetricContext.GetClusterEnergy (Cluster1);
					Try12=this->MetricContext.GetClusterEnergy (Cluster2);
					Try1=Try11+Try12;

					// Compute the energy when setting
					// I1 to the same cluster as I2;
					if ((*Size1 == 1)||(this->Self()->Derived::ThreadedConnexityConstraintProblem(I1,Edge,Val1,Val2,Process)==1))
						Try2 = 100000000.0;
					else
					{
						this->MetricContext.DeepCopy (Cluster1, Cluster21);
						this->MetricContext.DeepCopy (Cluster2, Cluster22);
						this->MetricContext.SubstractItemFromCluster (I1, Cluster21);
						this->MetricContext.AddItemToCluster (I1, Cluster22);
						this->MetricContext.ComputeClusterCentroid (Cluster21);
						this->MetricContext.ComputeClusterCentroid (Cluster22);
						this->MetricContext.ComputeClusterEnergy (Cluster21);
						this->MetricContext.ComputeClusterEnergy (Cluster22);

						Try21=this->MetricContext.GetClusterEnergy(Cluster21);
						Try22=this->MetricContext.GetClusterEnergy(Cluster22);
						Try2=Try21+Try22;
					}

					// Compute the energy when setting
					// I2 to the same cluster as I1;
					if ((*Size2 == 1) || (this->Self()->Derived::ThreadedConnexityConstraintProblem(I2, Edge, Val2, Val1, Process) == 1))
						Try3 = 1000000000.0;
					else
					{
						this->MetricContext.DeepCopy (Cluster1, Cluster31);
						this->MetricContext.DeepCopy (Cluster2, Cluster32);
						this->MetricContext.SubstractItemFromCluster (I2, Cluster32);
						this->MetricContext.AddItemToCluster (I2, Cluster31);
						this->MetricContext.ComputeClusterCentroid (Cluster31);
						this->MetricContext.ComputeClusterCentroid (Cluster32);
						this->MetricContext.ComputeClusterEnergy (Cluster31);
						this->MetricContext.ComputeClusterEnergy (Cluster32);

						Try31=this->MetricContext.GetClusterEnergy (Cluster31);
						Try32=this->MetricContext.GetClusterEnergy (Cluster32);
						Try3=Try31+Try32;
					}
					if ((Try1 <= Try2) && (Try1 <= Try3))
						Result = 1;
					else if ((Try2 < Try1) && (Try2 < Try3))
						Result = 2;
					else
						Result = 3;
				}
				else
				{
					Result = 1;
				}

				switch (Result)
				{
				case (1):
					// Don't do anything!
					this->vtkThreadedClustering::AddEdgeToProcess (Edge, Process);
					break;

				case (2):
					// Set I1 in the same cluster as I2
					this->Clustering->SetValue (I1, Val2);
					(*Size2)++;
					(*Size1)--;
					this->MetricContext.DeepCopy (Cluster21, Cluster1);
					this->MetricContext.DeepCopy (Cluster22, Cluster2);
					this->vtkThreadedClustering::AddItemRingToProcess (I1, Process, Process);
					data.NumberOfModifications++;
					this->ClustersLastModification[Val1] = this->NumberOfLoops;
					this->ClustersLastModification[Val2] = this->NumberOfLoops;
					break;

				case (3):
					// Set I2 in the same cluster as I1
					this->Clustering->SetValue (I2, Val1);
					(*Size1)++;
					(*Size2)--;
					this->MetricContext.DeepCopy (Cluster31, Cluster1);
					this->MetricContext.DeepCopy (Cluster32, Cluster2);
					this->vtkThreadedClustering::AddItemRingToProcess (I2, Process, Process);
					data.NumberOfModifications++;
					this->ClustersLastModification[Val1] = this->NumberOfLoops;
					this->ClustersLastModification[Val2] = this->NumberOfLoops;
				}
			}
		}
	}
	delete Cluster21;
	delete Cluster22;
	delete Cluster31;
	delete Cluster32;
}

template < class Metric, class Derived > void
	vtkThreadedClustering < Metric, Derived >::ExecuteProcessWithDistances (int Process)
{
	vtkIdType Edge, I1, I2;
	int Val1, Val2, *Size1, *Size2;
	ProcessData &data = this->ProcessInfos[ Process ];

	typename Metric::Cluster *Cluster1, *Cluster2;
	std::deque < int > *Queue;
	for (int CurrentQueue=0;CurrentQueue<this->PoolSize;CurrentQueue++)
	{
		Queue = &this->ProcessesPopQueues[Process][CurrentQueue];
		while (!Queue->empty())
		{
			// the queue is not empty : let's use it
			Edge = Queue->front ();
			Queue->pop_front ();
			this->Self()->Derived::GetEdgeItems (Edge, I1, I2);

			if ((this->EdgesLastLoop[Edge] !=this->RelativeNumberOfLoops) && (I2 >= 0))
			{
				this->EdgesLastLoop[Edge] = this->RelativeNumberOfLoops;
				Val1 = this->Clustering-> GetValue (I1);
				Val2 = this->Clustering-> GetValue (I2);
				if (Val2!=Val1)
				{

					data.NumberOfIterations++;
					if (Val1 == this->NumberOfClusters)
					{
						// I1 is not associated. Give it to the same cluster as I2
						this->MetricContext.AddItemToCluster(I1,&this->Clusters[Val2]);
						this->MetricContext.ComputeClusterCentroid(&this->Clusters[Val2]);
						this->MetricContext.ComputeClusterEnergy(&this->Clusters[Val2]);

						(*this->ClustersSizes->GetPointer (Val2))++;
						this->vtkThreadedClustering::AddItemRingToProcess (I1,Process, Process);
						this->Clustering->SetValue (I1,Val2);
						data.NumberOfModifications++;
					}
					else if (Val2 ==this->NumberOfClusters)
					{
						// I2 is not associated. Give it to the same cluster as I1
						this->MetricContext.AddItemToCluster(I2,&this->Clusters[Val1]);
						this->MetricContext.ComputeClusterCentroid(&this->Clusters[Val1]);
						this->MetricContext.ComputeClusterEnergy(&this->Clusters[Val1]);
						(*this->ClustersSizes->GetPointer (Val1))++;
						this->vtkThreadedClustering::AddItemRingToProcess (I2,Process,Process);
						this->Clustering->SetValue (I2,Val1);
						data.NumberOfModifications++;
					}
					else
					{
						int Result;

						// determine whether one of	the	two	adjacent clusters was modified,
						// or whether any of the clusters is freezed
						//	If not,	the	test is	useless, and the speed improved	:)
						if (((this->ClustersLastModification[Val1] >=this->NumberOfLoops-1)
						    || (this->ClustersLastModification[Val2]>=this->NumberOfLoops-1))
						    &&((this->IsClusterFreezed->GetValue(Val1)==0)
								&&(this->IsClusterFreezed->GetValue(Val2)==0)))
						{
							Cluster1=&this->Clusters[Val1];
							Cluster2=&this->Clusters[Val2];

							Size1=this->ClustersSizes->GetPointer(Val1);
							Size2=this->ClustersSizes->GetPointer(Val2);

							// Compute the initial energy
							double C1[3],C2[3];
							this->MetricContext.GetClusterCentroid(Cluster1,C1);
							this->MetricContext.GetClusterCentroid(Cluster2,C2);

							double P[3];
							// Compute the energy when setting I1 to the same cluster as I2;
							if ((*Size1!=1)&&(this->Self()->Derived::ThreadedConnexityConstraintProblem(I1,Edge,Val1,Val2,Process)==0))
							{
						//		this->MetricContext.GetItemCoordinates(I1,P);
								this->Self()->Derived::GetItemCoordinates(I1,P);
								if (vtkMath::Distance2BetweenPoints(P,C1)>
									vtkMath::Distance2BetweenPoints(P,C2))
								{
									Result=2;
									this->Clustering->SetValue(I1,Val2);
									(*Size2)++;
									(*Size1)--;
									this->MetricContext.AddItemToCluster(I1,Cluster2);
									this->MetricContext.SubstractItemFromCluster(I1,Cluster1);
									this->vtkThreadedClustering::AddItemRingToProcess (I1, Process, Process);
									data.NumberOfModifications++;
									this->ClustersLastModification[Val1]=this->NumberOfLoops;
									this->ClustersLastModification[Val2]=this->NumberOfLoops;
								}
							}
							if ((*Size2!=1)&&(this->Self()->Derived::ThreadedConnexityConstraintProblem(I2, Edge, Val2, Val1, Process)==0)
							&& (Result!=2))
							{
							// Compute the energy when setting I2 to the same cluster as I1;
							//	this->MetricContext.GetItemCoordinates(I2,P);
								this->Self()->Derived::GetItemCoordinates(I2,P);
								if (vtkMath::Distance2BetweenPoints(P,C1)<
									vtkMath::Distance2BetweenPoints(P,C2))
								{
									Result=3;
									this->Clustering->SetValue(I2,Val1);
									(*Size2)--;
									(*Size1)++;
									this->MetricContext.AddItemToCluster(I2,Cluster1);
									this->MetricContext.SubstractItemFromCluster(I2,Cluster2);
									this->vtkThreadedClustering::AddItemRingToProcess (I2, Process, Process);
									data.NumberOfModifications++;
									this->ClustersLastModification[Val1]=this->NumberOfLoops;
									this->ClustersLastModification[Val2]=this->NumberOfLoops;
								}
							}
						}
						else
						{
							Result = 1;
						}

						if (Result==1)
							this->vtkThreadedClustering::AddEdgeToProcess (Edge, Process);
					}
				}
			}
		}
	}
}


template < class Metric, class Derived >
int	vtkThreadedClustering < Metric, Derived >::ProcessOneLoop ()
{

	const int NumberOfRegions = this->PoolSize - 1;
	this->RegionClusters.resize (NumberOfRegions);
	vtkSMPTools::For (0, NumberOfRegions, 1, [this] (vtkIdType Begin, vtkIdType End)
	{
		for (vtkIdType Region = Begin; Region < End; Region++)
			this->CollectRegionClusters (static_cast<int> (Region));
	});
	this->ComputeRounds ();

	for ( auto &info : this->ProcessInfos ) info.NumberOfModifications=0;

	// the regions of a round share no state, so their order does not matter
	for (const std::vector<int> &Round : this->Rounds)
	{
		vtkSMPTools::For (0, static_cast<vtkIdType> (Round.size ()), 1,
			[this, &Round] (vtkIdType Begin, vtkIdType End)
		{
			for (vtkIdType i = Begin; i < End; i++)
				this->ExecuteRegion (Round[i]);
		});
	}

	// the edges between regions, alone
	ProcessData &data = this->ProcessInfos[this->PoolSize - 1];
	data.StartTime = vtkTimerLog::GetUniversalTime ();
	this->ExecuteProcess (this->PoolSize - 1);
	data.StopTime = vtkTimerLog::GetUniversalTime ();

	int NumberOfModifications=0;
	for ( auto &info : this->ProcessInfos )
		NumberOfModifications+=info.NumberOfModifications;

	this->DisplayThreadsTimings();
	return (NumberOfModifications);
}

template < class Metric, class Derived > void vtkThreadedClustering < Metric, Derived >::DisplayThreadsTimings()
{
	if (this->DisplayThreadsTimingsFlag)
	{
		ProcessData &firstData = this->ProcessInfos[0];
		double EarliestStart = firstData.StartTime;
		double LatestStop = firstData.StopTime;
		double MeasuredTime;
		double GlobalDuration;
		for (auto &info : this->ProcessInfos) {
			MeasuredTime = info.StartTime;
			if (MeasuredTime < EarliestStart)
				EarliestStart = MeasuredTime;
			MeasuredTime = info.StopTime;
			if (MeasuredTime > LatestStop)
				LatestStop = MeasuredTime;
		}
		GlobalDuration = LatestStop - EarliestStart;
		cout << "Processes duration:" << GlobalDuration <<" seconds"<<endl;

		cout << "Starts    (relative):";
		for (auto &info : this->ProcessInfos)
			cout << (int) (100.0 *(info.StartTime -EarliestStart) /GlobalDuration) << " ";

		cout << endl;
		cout << "Duration  (relative):";
		for (auto &info : this->ProcessInfos)
			cout << (int) (100.0 *(info.StopTime - info.StartTime)/GlobalDuration) << " ";

		cout << endl;
		cout << "Stops     (relative):";
		for (auto &info : this->ProcessInfos)
			cout << (int) (100.0 * (info.StopTime -	EarliestStart) / GlobalDuration) << " ";
		cout << endl;


		// Compute the clustering properties

		int MaxInstantNumberofIterations =
			firstData.NumberOfIterations - firstData.PreviousNumberOfIterations;
		int InstantNumber;
		int InstantNumberOfIeration = 0;
		for (auto &info : this->ProcessInfos) {
			InstantNumber =info.NumberOfIterations-info.PreviousNumberOfIterations;
			if (InstantNumber >MaxInstantNumberofIterations)
				MaxInstantNumberofIterations =InstantNumber;
			InstantNumberOfIeration +=info.NumberOfIterations-info.PreviousNumberOfIterations;
		}
		if (1)
		{
			cout << "Iterations(relative):";
			for (auto &info : this->ProcessInfos)
				cout << (int) (100.0 * ((double) (info.NumberOfIterations-info.PreviousNumberOfIterations)/((double)InstantNumberOfIeration)))<< " ";
			cout << endl;
			cout << "Iterations          :";
			for (auto &info : this->ProcessInfos)
				cout << info.NumberOfIterations<<" ";
			cout << endl;
			cout << "Modifications       :";
			for (auto &info : this->ProcessInfos)
				cout << info.NumberOfModifications<<" ";
			cout << endl;
		}

		for (auto &info : this->ProcessInfos)
			info.PreviousNumberOfIterations = info.NumberOfIterations;

	}
}

template < class Metric, class Derived > void vtkThreadedClustering < Metric, Derived >::SwapQueues ()
{
	std::deque < int >** Temp=this->ProcessesPopQueues;
	this->ProcessesPopQueues=this->ProcessesPushQueues;
	this->ProcessesPushQueues=Temp;
};


template < class Metric, class Derived > void vtkThreadedClustering < Metric, Derived >::Init ()
{
	vtkUniformClustering < Metric >::Init ();
	int i;

	// regions of at least ItemsPerRegion items and RegionClusterWidth clusters
	const vtkIdType NumberOfItems = this->GetNumberOfItems ();
	const vtkIdType ClusterSize = NumberOfItems / std::max (this->NumberOfClusters, 1);
	const vtkIdType RegionSize = std::max<vtkIdType> (
		{ this->ItemsPerRegion, RegionClusterWidth * ClusterSize, 1 });
	this->PoolSize = static_cast<int> (std::max<vtkIdType> (NumberOfItems / RegionSize, 1)) + 1;
	this->ComputeEdgesLayout ();

	// allocate statistics arrays
	this->ProcessInfos.resize(this->PoolSize);

	for ( auto &info : this->ProcessInfos ) {
		info.NumberOfIterations = 0;
		info.PreviousNumberOfIterations = 0;
		info.NumberOfModifications = 0;
	}

	// allocate queues
	this->ProcessesQueues1=new std::deque<int>*[this->PoolSize];
	this->ProcessesQueues2=new std::deque<int>*[this->PoolSize];
	for (i = 0; i < this->PoolSize; i++)
	{
		this->ProcessesQueues1[i]=new std::deque<int>[this->PoolSize];
		this->ProcessesQueues2[i]=new std::deque<int>[this->PoolSize];
	}

	// set initial queues state
	this->ProcessesPushQueues=this->ProcessesQueues1;
	this->ProcessesPopQueues=this->ProcessesQueues2;

}

template < class Metric, class Derived > void vtkThreadedClustering < Metric, Derived >::ComputeEdgesLayout ()
{
	vtkIntArray *TempClustering=0;
	vtkTimerLog *Timer=vtkTimerLog::New();
	Timer->StartTimer();
	switch (this->EdgesLayoutComputingType)
	{
	case 1:
		if (this->ConsoleOutput)
			cout<<"Creating "<<this->PoolSize-1<<" Regions for thread pooling...";
		TempClustering = vtkIntArray::New ();
		TempClustering->SetNumberOfValues (this->GetNumberOfItems ());
		this->ComputeInitialRandomSampling (0, TempClustering,this->PoolSize-1);
		this->FillHolesInClustering (TempClustering);
		if (this->ConsoleOutput)
			cout<<"done"<<endl;
		break;

	case 2:
	{
		/*
		 * if (this->ClusteringType==0)
		 * {
		 * vtkTrianglesUniformClustering *Clustering=vtkTrianglesUniformClustering::New();
		 * Clustering->SetInput(this->Input);
		 * Clustering->SetDisplay(0);
		 * Clustering->Setthis->NumberOfClusters(this->NumberOfProcesses);
		 * TempClustering=Clustering->ProcessClustering();
		 * }
		 * else
		 * {
		 * vtkVerticesUniformClustering *Clustering=vtkVerticesUniformClustering::New();
		 * Clustering->SetInput(this->Input);
		 * Clustering->SetDisplay(0);
		 * Clustering->SetNumberOfClusters(this->NumberOfProcesses);
		 * TempClustering=Clustering->ProcessClustering();
		 * }
		 * break;
		 */
	}

	case 3:
		TempClustering = vtkIntArray::New ();
		TempClustering->SetNumberOfValues (this->GetNumberOfItems ());
		this->ComputeMeshSlices (this->PoolSize-1, TempClustering);
	}

	int C1;
	int I1;
	int NumberOfProblems = 0;
	int i, j;
	int *EdgesLayout = new int[this->GetNumberOfEdges ()];

	this->RegionNeighbours.assign (this->PoolSize - 1, std::vector<int> ());
	vtkIdList *VList = vtkIdList::New ();
	for (i = 0; i < this->GetNumberOfEdges (); i++)
	{
		this->GetEdgeItemsSure (i, VList);
		I1 = VList->GetId (0);
		C1 = TempClustering->GetValue (I1);
		if ((C1 >= this->PoolSize) || (C1 < 0))
		{
			NumberOfProblems++;
			C1 = this->PoolSize - 1;
		}
		for (j = 0; j < VList->GetNumberOfIds (); j++)
		{
			if (TempClustering->GetValue (VList->GetId (j)) != C1)
			{
				C1 = this->PoolSize - 1;
			}
		}
		EdgesLayout[i] = C1;

		// the regions an edge between regions joins are adjacent
		if (C1 == this->PoolSize - 1)
		{
			for (j = 0; j < VList->GetNumberOfIds (); j++)
			{
				int R1 = TempClustering->GetValue (VList->GetId (j));
				if ((R1 < 0) || (R1 >= this->PoolSize - 1)) continue;
				for (int k = j + 1; k < VList->GetNumberOfIds (); k++)
				{
					int R2 = TempClustering->GetValue (VList->GetId (k));
					if ((R2 < 0) || (R2 >= this->PoolSize - 1) || (R2 == R1)) continue;
					this->RegionNeighbours[R1].push_back (R2);
					this->RegionNeighbours[R2].push_back (R1);
				}
			}
		}
	}
	for (std::vector<int> &Neighbours : this->RegionNeighbours)
	{
		std::sort (Neighbours.begin (), Neighbours.end ());
		Neighbours.erase (std::unique (Neighbours.begin (), Neighbours.end ()), Neighbours.end ());
	}
	if (NumberOfProblems != 0)
		cout << "Partitionning for multithreading: " << NumberOfProblems <<
			" Items with no good association..." << endl;
	VList->Delete ();
	this->EdgesProcess = std::vector<int>( EdgesLayout, EdgesLayout+this->GetNumberOfEdges());
	Timer->StopTimer();
	cout<<"Multithreading layout computed in :"<<Timer->GetElapsedTime()<<" seconds"<<endl;
	Timer->Delete();
	TempClustering->Delete();
}


template < class Metric, class Derived >
void vtkThreadedClustering < Metric, Derived >::ComputeMeshSlices (int NumberOfSlices,vtkIntArray * Clustering)
{
	int i,j;
	int Granularity = 10000;
	double *Weights = new double[Granularity];
	for (i = 0; i < Granularity; i++)
		Weights[i] = 0;

	double P[3];

	int *Slices = new int[Granularity];

	// compute the bounding box to choose the best direction
	double Bounds[6];
	this->GetItemCoordinates (0, P);
	Bounds[0]=P[0];
	Bounds[1]=P[0];
	Bounds[2]=P[1];
	Bounds[3]=P[1];
	Bounds[4]=P[2];
	Bounds[5]=P[2];

	for (i = 0; i < this->GetNumberOfItems (); i++)
	{
		this->GetItemCoordinates (i, P);
		for (j=0;j<3;j++)
		{
			if (Bounds[2*j]>P[j])
				Bounds[2*j]=P[j];
			if (Bounds[2*j+1]<P[j])
				Bounds[2*j+1]=P[j];
		}
	}

	double DeltaMax = 0;
	double Delta;
	int Direction = -1;
	double CoordMax=-1, CoordMin=100000000;
	for (i = 0; i < 3; i++)
	{
		Delta = fabs (Bounds[2 * i + 1] - Bounds[2 * i]);
		if (DeltaMax < Delta)
		{
			DeltaMax = Delta;
			Direction = i;
			CoordMax = Bounds[2 * i + 1];
			CoordMin = Bounds[2 * i];
		}
	}

	double Step = (CoordMax - CoordMin) / (double) (Granularity - 1);
	int Index;
	for (i = 0; i < this->GetNumberOfItems (); i++)
	{
		this->GetItemCoordinates (i, P);
		Index = (int) floor ((P[Direction] - CoordMin) / Step);
		if ((Index < 0) || (Index >= Granularity))
		{
			cout << "Item " << i << " out of bounds ! Index=" <<
				Index << endl;
			cout << "Contact developpers...." << endl;
			Index = 0;
		}
		Weights[Index] += this->MetricContext.GetItemWeight(i);
		Clustering->SetValue (i, Index);
	}

	// Compute the average weight for each slice
	double AverageWeight = 0;
	for (i = 0; i < Granularity; i++)
		AverageWeight += Weights[i];
	AverageWeight = AverageWeight / (double) NumberOfSlices;

	// assign the slices
	double SWeight = 0;
	Index = 0;
	for (i = 0; i < Granularity; i++)
	{
		Slices[i] = Index;
		SWeight += Weights[i];
		if (SWeight > AverageWeight)
		{
			SWeight -= AverageWeight;
			Index++;
		}
	}

	// update the clustering
	for (i = 0; i < this->GetNumberOfItems (); i++)
	{
		Clustering->SetValue (i, Slices[Clustering->GetValue (i)]);
	}

	delete[]Slices;
	delete[]Weights;
}


template < class Metric, class Derived >
	vtkThreadedClustering < Metric, Derived >::vtkThreadedClustering ()
{
	this->EdgesLayoutComputingType = 1;
	this->DisplayThreadsTimingsFlag = 0;

	this->Timer = vtkTimerLog::New ();

	this->ItemsPerRegion = 512;

	this->ProcessesQueues1=0;
	this->ProcessesQueues2=0;

	this->NumberOfClusters=0;
}

template < class Metric, class Derived >
	vtkThreadedClustering < Metric, Derived >::~vtkThreadedClustering ()
{
	this->Timer->Delete ();

	if (this->ProcessesQueues1)
	{
		for (int i=0;i<PoolSize;i++)
		{
			delete [] this->ProcessesQueues1[i];
			delete [] this->ProcessesQueues2[i];
		}
		delete [] this->ProcessesQueues1;
		delete [] this->ProcessesQueues2;

	}
}

#endif
