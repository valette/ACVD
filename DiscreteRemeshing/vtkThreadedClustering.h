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

#ifndef _VTKTHREADEDCLUSTERING_H_
#define _VTKTHREADEDCLUSTERING_H_

#include <vtkMath.h>
#include <vtkTimerLog.h>
#include <vtkMultiThreader.h>
#include <vtkPriorityQueue.h>
#include <mutex>

#include "vtkUniformClustering.h"

// Class derived from vtkUniformClustering


// this switch enables/disable the use of mutexes to ensure thread safety
// although the risk of collision is very low when processing large meshes with a low number of threads,
// this switch should be ALWAYS ON
#define THREADSAFECLUSTERING

template < class Metric > class vtkThreadedClustering:public vtkUniformClustering <	Metric >
{
public:

	/// Sets the number of desired threads for multithreading.
	/// by default, the number of threads is the number of CPUS
	void SetNumberOfThreads (int N)
	{
		this->NumberOfThreads=N;
	};

	/// returns the number of threads used for the clustering.
	int GetNumberOfThreads ()
	{
		return (this->NumberOfThreads);
	};
	
	/// sets the pooling ratio (the number of pool jobs will be equal to NumberOfThreads*PoolingRatio+1
	/// default value : 5
	void SetPoolingRatio (int Ratio)
	{
		this->PoolingRatio=Ratio;
	};
	
	virtual vtkIntArray *ProcessClustering(vtkIdList *List=0)
	{
		this->vtkUniformClustering<Metric>::ProcessClustering();
		int NumberOfLockingCollisions=0;
		for (int i=0;i<this->NumberOfThreads+1;i++)
			NumberOfLockingCollisions+=this->threadInfos[ i ].NumberOfLockingCollisions;

		if (NumberOfLockingCollisions!=0)
		cout<<NumberOfLockingCollisions<<" locking collisions detected"<<endl;
			else
		cout<<"No locking collision detected"<<endl;
		return (this->Clustering);
	};

protected:

	// the constructor
	vtkThreadedClustering ();
	
	// the destructor
	~vtkThreadedClustering ();

	// Executes one process from the Thread pool
	void ExecuteProcess (int Process, int Thread);
	void ExecuteProcessWithDistances (int Process,int Thread);

	// The static function used for clustering
	static VTK_THREAD_RETURN_TYPE MyMainForClustering (void *arg);

	// virtual function. Might be implemented in derived classes for speed issues
	virtual void AddItemRingToProcess (vtkIdType Item, int ProcessId, int Thread)
	{
		vtkIdList *EList=this->threadInfos[Thread].itemList;
		this->GetItemEdges(Item,EList);
		for (int i=0;i<EList->GetNumberOfIds();i++)
			this->AddEdgeToProcess(EList->GetId(i),ProcessId);
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
					  EdgesProcess[Edge]][ProcessId].push (Edge);
	};

	/// The number of threads
	int NumberOfThreads;

	/// The size of the pool
	int PoolSize;
	
	// Parameter defining the number of thread jobs : NThreadJobs=NumberOfThreads*PoolingRatio+1
	int PoolingRatio;

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

	class ThreadData {
		public:
		int PreviousNumberOfIterations;
		int NumberOfIterations;
		int NumberOfModifications;
		int NumberOfLockingCollisions;
		double StartTime;
		double StopTime;
		vtkIdList *itemList;
		ThreadData() { itemList = vtkIdList::New(); };
		~ThreadData() { itemList->Delete(); };
	};

	std::vector< ThreadData > threadInfos;

#ifdef THREADSAFECLUSTERING
	// one mutex for each cluster to ensure thread-safety
	std::vector < std::mutex > ClustersLocks;
#endif

	// Context used to allocate processes in the pool
	// **********************************************
	std::mutex PoolAllocationLock,PoolAllocationLock2;
	vtkPriorityQueue *PoolQueue1,*PoolQueue2;

	// The push and pop queues used to track the edges laying between different clusters
	std::queue <int>**ProcessesQueues1;
	std::queue <int>**ProcessesQueues2;

	std::queue <int>**ProcessesPushQueues;
	std::queue <int>**ProcessesPopQueues;
	
	// this method swaps the pop queues with the push queues.
	void SwapQueues ();
};

template < class Metric > void vtkThreadedClustering < Metric >::FillQueuesFromClustering ()
{
	int i, j;
	vtkIdType I1, I2;
	std::queue < int >*PQueue;
	vtkIdType Edge;

	// empty the pushing queues
	for (i = 0; i < this->PoolSize; i++)
	{
		for (j = 0; j < this->PoolSize; j++)
		{
			PQueue = &this->ProcessesPushQueues[i][j];
			while (PQueue->size () != 0)
			{
				Edge = PQueue->front ();
				PQueue->pop ();
			}
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
					[this->EdgesProcess[i]].push (i);
			}
		}
	}
	
	// fill the thread pool queue
	for (i=0;i<this->PoolSize-1;i++)
		this->PoolQueue1->Insert(i,i);	
}


template < class Metric >	VTK_THREAD_RETURN_TYPE vtkThreadedClustering 
<Metric >::MyMainForClustering (void *arg)
{
	vtkMultiThreader::ThreadInfo * Info =
		(vtkMultiThreader::ThreadInfo *) arg;

	vtkThreadedClustering < Metric > *Clustering =
		(vtkThreadedClustering < Metric > *)Info->UserData;

	int MyId = Info->ThreadID;
	ThreadData &data = Clustering->threadInfos[ MyId ];

	data.StartTime = Clustering->Timer->GetUniversalTime ();

	while (1)
	{
		Clustering->PoolAllocationLock.lock();
		int Process=Clustering->PoolQueue1->Pop();
		Clustering->PoolAllocationLock.unlock();
				
		if ((Process<0))
			break;

		double Time;

		Time = Clustering->Timer->GetUniversalTime ();
		if (Clustering->MinimizeUsingEnergy)
			Clustering->ExecuteProcess (Process,MyId);
		else
			Clustering->ExecuteProcessWithDistances (Process,MyId);
		Time = Clustering->Timer->GetUniversalTime ()-Time;
		Clustering->PoolAllocationLock2.lock();
		Clustering->PoolQueue2->Insert (-Time,Process);
		Clustering->PoolAllocationLock2.unlock();
	}

	data.StopTime = Clustering->Timer->GetUniversalTime ();
	return VTK_THREAD_RETURN_VALUE;
}

template < class Metric > void
	vtkThreadedClustering < Metric >::ExecuteProcess (int Process,int Thread)
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
	ThreadData &data = this->threadInfos[ Thread ];

	std::queue < int > *Queue;
	for (int CurrentQueue=0;CurrentQueue<this->PoolSize;CurrentQueue++)
	{
		Queue = &this->ProcessesPopQueues[Process][CurrentQueue];
		while (!Queue->empty())
		{
			// the queue is not empty : let's use it
			Edge = Queue->front ();
			Queue->pop ();
			this->GetEdgeItems (Edge, I1, I2);

			if ((this->EdgesLastLoop[Edge] ==this->RelativeNumberOfLoops) || (I2 < 0)) continue;
			this->EdgesLastLoop[Edge] = this->RelativeNumberOfLoops;
			Val1 = this->Clustering-> GetValue (I1);
			Val2 = this->Clustering-> GetValue (I2);
			if (Val2==Val1) continue;
#ifdef THREADSAFECLUSTERING
			// get the lock on the clusters
			if (Val1<Val2)
			{
				if (!this->ClustersLocks[Val1].try_lock())
				{
					data.NumberOfLockingCollisions++;
					this->ClustersLocks[Val1].lock();
				}
				if (!this->ClustersLocks[Val2].try_lock())
				{
					data.NumberOfLockingCollisions++;
					this->ClustersLocks[Val2].lock();
				}
			}
			else
			{
				if (!this->ClustersLocks[Val2].try_lock())
				{
					data.NumberOfLockingCollisions++;
					this->ClustersLocks[Val2].lock();
				}
				if (!this->ClustersLocks[Val1].try_lock())
				{
					data.NumberOfLockingCollisions++;
					this->ClustersLocks[Val1].lock();
				}
			}

#endif

			data.NumberOfIterations++;
			if (Val1 == this->NumberOfClusters)
			{
				// I1 is not associated. Give it to the same cluster as I2
				this->MetricContext.AddItemToCluster(I1,&this->Clusters[Val2]);
				this->MetricContext.ComputeClusterCentroid(&this->Clusters[Val2]);
				this->MetricContext.ComputeClusterEnergy(&this->Clusters[Val2]);
				(*this->ClustersSizes->GetPointer (Val2))++;
				this->AddItemRingToProcess (I1,Process, Thread);
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
				this->AddItemRingToProcess (I2,Process,Thread);
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
					if ((*Size1 == 1)||(this->ThreadedConnexityConstraintProblem(I1,Edge,Val1,Val2,Thread)==1))
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
					if ((*Size2 == 1) || (this->ThreadedConnexityConstraintProblem (I2, Edge, Val2, Val1, Thread) == 1))
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
					this->AddEdgeToProcess (Edge, Process);
					break;

				case (2):
					// Set I1 in the same cluster as I2
					this->Clustering->SetValue (I1, Val2);
					(*Size2)++;
					(*Size1)--;
					this->MetricContext.DeepCopy (Cluster21, Cluster1);
					this->MetricContext.DeepCopy (Cluster22, Cluster2);
					this->AddItemRingToProcess (I1, Process, Thread);
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
					this->AddItemRingToProcess (I2, Process, Thread);
					data.NumberOfModifications++;
					this->ClustersLastModification[Val1] = this->NumberOfLoops;
					this->ClustersLastModification[Val2] = this->NumberOfLoops;
				}
			}
#ifdef THREADSAFECLUSTERING
		// release the lock on the clusters
		this->ClustersLocks[Val1].unlock();
		this->ClustersLocks[Val2].unlock();
#endif	
		}
	}
	delete Cluster21;
	delete Cluster22;
	delete Cluster31;
	delete Cluster32;	
}

template < class Metric > void
	vtkThreadedClustering < Metric >::ExecuteProcessWithDistances (int Process,int Thread)
{
	vtkIdType Edge, I1, I2;
	int Val1, Val2, *Size1, *Size2;
	ThreadData &data = this->threadInfos[ Thread ];

	typename Metric::Cluster *Cluster1, *Cluster2;
	std::queue < int > *Queue;
	for (int CurrentQueue=0;CurrentQueue<this->PoolSize;CurrentQueue++)
	{
		Queue = &this->ProcessesPopQueues[Process][CurrentQueue];
		while (!Queue->empty())
		{
			// the queue is not empty : let's use it
			Edge = Queue->front ();
			Queue->pop ();
			this->GetEdgeItems (Edge, I1, I2);

			if ((this->EdgesLastLoop[Edge] !=this->RelativeNumberOfLoops) && (I2 >= 0))
			{
				this->EdgesLastLoop[Edge] = this->RelativeNumberOfLoops;
				Val1 = this->Clustering-> GetValue (I1);
				Val2 = this->Clustering-> GetValue (I2);
				if (Val2!=Val1)
				{
#ifdef THREADSAFECLUSTERING
					// get the lock on the clusters
					if (Val1<Val2)
					{
						if (!this->ClustersLocks[Val1].try_lock())
						{
							data.NumberOfLockingCollisions++;
							this->ClustersLocks[Val1].lock();
						}
						if (!this->ClustersLocks[Val2].try_lock())
						{
							data.NumberOfLockingCollisions++;
							this->ClustersLocks[Val2].lock();
						}
					}
					else
					{
						if (!this->ClustersLocks[Val2].try_lock())
						{
							data.NumberOfLockingCollisions++;
							this->ClustersLocks[Val2].lock();
						}
						if (!this->ClustersLocks[Val1].try_lock())
						{
							data.NumberOfLockingCollisions++;
							this->ClustersLocks[Val1].lock();
						}
					}

#endif

					data.NumberOfIterations++;
					if (Val1 == this->NumberOfClusters)
					{
						// I1 is not associated. Give it to the same cluster as I2
						this->MetricContext.AddItemToCluster(I1,&this->Clusters[Val2]);
						this->MetricContext.ComputeClusterCentroid(&this->Clusters[Val2]);
						this->MetricContext.ComputeClusterEnergy(&this->Clusters[Val2]);

						(*this->ClustersSizes->GetPointer (Val2))++;
						this->AddItemRingToProcess (I1,Process, Thread);
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
						this->AddItemRingToProcess (I2,Process,Thread);
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
							if ((*Size1!=1)&&(this->ThreadedConnexityConstraintProblem(I1,Edge,Val1,Val2,Thread)==0))
							{
						//		this->MetricContext.GetItemCoordinates(I1,P);
								this->GetItemCoordinates(I1,P);
								if (vtkMath::Distance2BetweenPoints(P,C1)>
									vtkMath::Distance2BetweenPoints(P,C2))
								{
									Result=2;
									this->Clustering->SetValue(I1,Val2);
									(*Size2)++;
									(*Size1)--;
									this->MetricContext.AddItemToCluster(I1,Cluster2);
									this->MetricContext.SubstractItemFromCluster(I1,Cluster1);
									this->AddItemRingToProcess (I1, Process, Thread);
									data.NumberOfModifications++;
									this->ClustersLastModification[Val1]=this->NumberOfLoops;
									this->ClustersLastModification[Val2]=this->NumberOfLoops;
								}
							}
							if ((*Size2!=1)&&(this->ThreadedConnexityConstraintProblem (I2, Edge, Val2, Val1, Thread)==0)
							&& (Result!=2))
							{
							// Compute the energy when setting I2 to the same cluster as I1;
							//	this->MetricContext.GetItemCoordinates(I2,P);
								this->GetItemCoordinates(I2,P);
								if (vtkMath::Distance2BetweenPoints(P,C1)<
									vtkMath::Distance2BetweenPoints(P,C2))
								{
									Result=3;
									this->Clustering->SetValue(I2,Val1);
									(*Size2)--;
									(*Size1)++;
									this->MetricContext.AddItemToCluster(I2,Cluster1);
									this->MetricContext.SubstractItemFromCluster(I2,Cluster2);
									this->AddItemRingToProcess (I2, Process, Thread);
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
							this->AddEdgeToProcess (Edge, Process);
					}
#ifdef THREADSAFECLUSTERING
				// release the lock on the clusters
				this->ClustersLocks[Val1].unlock();
				this->ClustersLocks[Val2].unlock();
#endif
				}
			}
		}
	}
}


template < class Metric >
int	vtkThreadedClustering < Metric >::ProcessOneLoop ()
{

#ifdef THREADSAFECLUSTERING
	if ( this->ClustersLocks.size() != ( this->NumberOfClusters+1 ) )
		this->ClustersLocks=std::vector< std::mutex >(this->NumberOfClusters+1);
#endif

	vtkMultiThreader *Threader=vtkMultiThreader::New();
	Threader->SetSingleMethod (MyMainForClustering, (void *) this);
	Threader->SetNumberOfThreads (this->NumberOfThreads);
		
	for ( auto &info : this->threadInfos ) info.NumberOfModifications=0;
	Threader->SingleMethodExecute ();

	this->threadInfos[this->NumberOfThreads].StartTime = this->Timer->GetUniversalTime ();
	this->ExecuteProcess (this->PoolSize - 1,this->NumberOfThreads);
	this->threadInfos[this->NumberOfThreads].StopTime = this->Timer->GetUniversalTime ();

	vtkPriorityQueue *Queue=this->PoolQueue1;
	this->PoolQueue1=this->PoolQueue2;
	this->PoolQueue2=Queue;
	Threader->Delete();

	int NumberOfModifications=0;
	for ( auto &info : this->threadInfos )
		NumberOfModifications+=info.NumberOfModifications;

	this->DisplayThreadsTimings();
	return (NumberOfModifications);
}

template < class Metric > void vtkThreadedClustering < Metric >::DisplayThreadsTimings()
{
	if (this->DisplayThreadsTimingsFlag)
	{
		ThreadData &firstData = this->threadInfos[0];
		double EarliestStart = firstData.StartTime;
		double LatestStop = firstData.StopTime;
		double MeasuredTime;
		double GlobalDuration;
		for (auto &info : this->threadInfos) {
			MeasuredTime = info.StartTime;
			if (MeasuredTime < EarliestStart)
				EarliestStart = MeasuredTime;
			MeasuredTime = info.StopTime;
			if (MeasuredTime > LatestStop)
				LatestStop = MeasuredTime;
		}
		GlobalDuration = LatestStop - EarliestStart;
		cout << "Threads duration:" << GlobalDuration <<" seconds"<<endl;

		cout << "Starts    (relative):";
		for (auto &info : this->threadInfos)
			cout << (int) (100.0 *(info.StartTime -EarliestStart) /GlobalDuration) << " ";

		cout << endl;
		cout << "Duration  (relative):";
		for (auto &info : this->threadInfos)
			cout << (int) (100.0 *(info.StopTime - info.StartTime)/GlobalDuration) << " ";

		cout << endl;
		cout << "Stops     (relative):";
		for (auto &info : this->threadInfos)
			cout << (int) (100.0 * (info.StopTime -	EarliestStart) / GlobalDuration) << " ";
		cout << endl;


		// Compute the clustering properties

		int MaxInstantNumberofIterations =
			firstData.NumberOfIterations - firstData.PreviousNumberOfIterations;
		int InstantNumber;
		int InstantNumberOfIeration = 0;
		for (auto &info : this->threadInfos) {
			InstantNumber =info.NumberOfIterations-info.PreviousNumberOfIterations;
			if (InstantNumber >MaxInstantNumberofIterations)
				MaxInstantNumberofIterations =InstantNumber;
			InstantNumberOfIeration +=info.NumberOfIterations-info.PreviousNumberOfIterations;
		}
		if (1)
		{
			cout << "Iterations(relative):";
			for (auto &info : this->threadInfos)
				cout << (int) (100.0 * ((double) (info.NumberOfIterations-info.PreviousNumberOfIterations)/((double)InstantNumberOfIeration)))<< " ";
			cout << endl;
			cout << "Iterations          :";
			for (auto &info : this->threadInfos)
				cout << info.NumberOfIterations<<" ";
			cout << endl;			
			cout << "Modifications       :";
			for (auto &info : this->threadInfos)
				cout << info.NumberOfModifications<<" ";
			cout << endl;			
		}

		for (auto &info : this->threadInfos)
			info.PreviousNumberOfIterations = info.NumberOfIterations;

	}
}

template < class Metric > void vtkThreadedClustering < Metric >::SwapQueues ()
{
	std::queue < int >** Temp=this->ProcessesPopQueues;
	this->ProcessesPopQueues=this->ProcessesPushQueues;
	this->ProcessesPushQueues=Temp;
};


template < class Metric > void vtkThreadedClustering < Metric >::Init ()
{
	vtkUniformClustering < Metric >::Init ();
	int i;
	
	this->PoolSize=this->PoolingRatio*NumberOfThreads+1;
	this->ComputeEdgesLayout ();

	// allocate statistics arrays
	this->threadInfos.resize(this->NumberOfThreads+1);

	for ( auto &info : this->threadInfos ) {
		info.NumberOfIterations = 0;
		info.PreviousNumberOfIterations = 0;
		info.NumberOfModifications = 0;
		info.NumberOfLockingCollisions = 0;
	}

	// allocate queues
	this->ProcessesQueues1=new std::queue<int>*[this->PoolSize];
	this->ProcessesQueues2=new std::queue<int>*[this->PoolSize];
	for (i = 0; i < this->PoolSize; i++)
	{
		this->ProcessesQueues1[i]=new std::queue<int>[this->PoolSize];
		this->ProcessesQueues2[i]=new std::queue<int>[this->PoolSize];		
	}
	
	// set initial queues state
	this->ProcessesPushQueues=this->ProcessesQueues1;
	this->ProcessesPopQueues=this->ProcessesQueues2;

}

template <class Metric> void vtkThreadedClustering<Metric>::ComputeEdgesLayout ()
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

	
template <class Metric> 
void vtkThreadedClustering<Metric>::ComputeMeshSlices (int NumberOfSlices,vtkIntArray * Clustering)
{
	int i,j;
	int Granularity = 10000;
	double *Weights = new double[Granularity];
	for (i = 0; i < Granularity; i++)
		Weights[i] = 0;

	double P[3];

	int *Slices = new int[Granularity];

	// compute the bouding box to choose the best direction
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


template < class Metric >
	vtkThreadedClustering < Metric >::vtkThreadedClustering ()
{
	this->EdgesLayoutComputingType = 1;
	this->DisplayThreadsTimingsFlag = 0;

	this->Timer = vtkTimerLog::New ();

	// autoset the number of threads depending on the number of processors
	vtkMultiThreader *Threader=vtkMultiThreader::New();
	this->SetNumberOfThreads(Threader->GetNumberOfThreads());
	Threader->Delete();
	
	this->PoolingRatio=5;
	
	this->PoolQueue1=vtkPriorityQueue::New();
	this->PoolQueue2=vtkPriorityQueue::New();
	
	this->ProcessesQueues1=0;
	this->ProcessesQueues2=0;
	
	this->NumberOfClusters=0;
}

template < class Metric >
	vtkThreadedClustering < Metric >::~vtkThreadedClustering ()
{
	this->Timer->Delete ();
	this->PoolQueue1->Delete();
	this->PoolQueue2->Delete();

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
