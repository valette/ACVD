/*=========================================================================

  Program:   Uniform Clustering (abstract class)
  Module:    vtkUniformClustering.h
  Language:  C++
  Date:      2004/09
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

#ifndef _VTKUNIFORMCLUSTERING_H_
#define _VTKUNIFORMCLUSTERING_H_

#include <algorithm>
#include <vtkCommand.h>
#include <vtkMath.h>
#include <vtkCellData.h>
#include <vtkTimerLog.h>
#include "RenderWindow.h"

/// A Class to process uniform clustering, Implemented from the paper:
/// "Approximated Centroidal Voronoi Diagrams for Uniform Polygonal Mesh Coarsening"
/// [Valette & Chassery, Eurographics 2004]
/// NOTE : this is a pure abstract class, only derived class can be used.
/// See vtkSurfaceClustering for a working example

template <class Metric, typename EdgeType=vtkIdType> class vtkUniformClustering:
public vtkObject
{
public:
	
	/// Sets the number of clusters to create
	void SetNumberOfClusters(int N){this->NumberOfClusters=N;};
	
	/// Main class call: processes the clustering. List is the list of items to cluster. If you want to process
	/// all the items (which is the case in 99.99% of the cases), leave this parameter empty
	vtkIntArray *ProcessClustering(vtkIdList *List=0);
	
	/// Returns the Clustering
	vtkIntArray *GetClustering();
	
	// Sets On/Off the fast initialization by unconstrained clustering (default : Off)
	void SetConstrainedInitialization(int C) {this->UnconstrainedInitializationFlag=C;};
			
	/// Enables/Disables graphical display
	/// 0: no display  1: display after convergence   2: display during convergence (might be slow)
	/// 3: display during clustering and capture images
	/// default:1
	void SetDisplay(int d){this->Display=d;};
	
	/// Sets the Anchor Renderwindow of the class
	/// this is usefull when several windows have to be linked outside this class
	void SetAnchorRenderWindow(RenderWindow *Window) {this->AnchorRenderWindow=Window;};
	
	/// Enables/Disables console output (text mode)  while processing
	/// 0: off  1:on  Default:0
	void SetConsoleOutput(int d){this->ConsoleOutput=d;};
	
	/// Sets the maximun number of convergences while minimizing the energy term
	/// can speed up things a bit (usually not needed)
	void SetMaxNumberOfConvergences(int n){this->MaxNumberOfConvergences=n;};
	
	/// Sets the maximum number of loops during minimzation
	/// can be usefull when there are always clusters with disconnected components (VERY VERY RARE)
	void SetMaxNumberOfLoops(int n){this->MaxNumberOfLoops=n;};
		
	/// Determines how the initial sampling is done
	/// 0: each cluster gets one randomly picked item
	/// 1: the clusters are initialized according to the density function
	void SetInitialSamplingType (int Type) {this->InitialSamplingType=Type;};

	/// Defines an initial clustering 
	void SetInitialClustering(vtkIntArray *Clustering);
	
	/// Returns the Number of clusters
	int GetNumberOfClusters() {return this->NumberOfClusters;};

	/// Returns a pointer to the metric. Usefull to Set metric-dependent options (Gradation, etc...)
	Metric* GetMetric() {return (&this->MetricContext);};

	/// Returns a pointer to a given cluster
	typename Metric::Cluster *GetCluster(int C) {return (this->Clusters+C);};
	
	/// Sets On/Off the writing of the energy evolution during the clustering to a file called "energy.txt"
	void SetWriteToGlobalEnergyLog(int S) {this->ComputeAndSaveEnergy=S;};
	
	/// Defines the output directory where thr log file will be written
	void SetOutputDirectory(char *S) {this->OutputDirectory=S;};

	/// Pure virtual method to be implemented in derived classes
	virtual int GetNumberOfItems()=0;

protected:

	/// virtual function derived in the threaded version.
	virtual void SwapQueues () {};

	/// Virtual method that can be overcharged in derived classes when speed is an issue
	virtual void AddItemRingToProcess(vtkIdType Item);

	/// Pure virtual method to be implemented in derived classes
	virtual void GetItemNeighbours(vtkIdType Item,vtkIdList *IList)=0;

	/// Pure virtual method to be implemented in derived classes
	virtual EdgeType GetNumberOfEdges()=0;

	/// Pure virtual method to be implemented in derived classes
	virtual void GetEdgeItems(EdgeType Edge, vtkIdType &I1,vtkIdType &I2)=0;

	/// Pure virtual method to be implemented in derived classes
	virtual void GetItemEdges(vtkIdType Item, vtkIdList *EList)=0;

	/// Pure virtual method to be implemented in derived classes
	virtual int ConnexityConstraintProblem(vtkIdType Item,EdgeType Edge,vtkIdType Cluster, vtkIdType Cluster2)=0;

	// Virtual method used to add constraints when cleaning the clustering
	virtual bool IsClusterCleanable(vtkIdType Cluster)
		{return (true);};

	/// Pure virtual method to be implemented in derived classes
	virtual void GetItemCoordinates(vtkIdType Item,double *P)=0;

	/// displays the clustering in the RenderWindow, possibly captures the images depending on Display options
	virtual void Snapshot()=0;

	/// method to create the render window(s)
	virtual void CreateWindows()=0;
	
	/// method to enable window interaction
	virtual void InteractWithClusteringWindow(){};

	/// Method to build the metric (pure virtual function)
	virtual void BuildMetric()=0;
	
	/// Method mostly usefull when clustering volumes with boundaries constraints : 
	/// this method returns whether an item is on the boundary between objects
	/// the method should return -1 when the item is not on a boundary.
	virtual int GetItemType(vtkIdType Item){return (1);};

	/// The metric itself, containing all the clusters-items interactions
	Metric MetricContext;

	/// The clusters, metric-dependent
	typename Metric::Cluster *Clusters;

	/// the Anchor render window used to link all the windows
	RenderWindow *AnchorRenderWindow;

	/// the constructor
	vtkUniformClustering(); 

	/// the destructor
	~vtkUniformClustering(); 
	
	/// Parameter indicating if the connexity constraint for the clusters has to be respected
	int ConnexityConstraint;

	/// The Clustrering (one value for each Item)
	vtkIntArray *Clustering;
	
	/// Possibly points to a user-defined initial clustering when one was defined by SetInitialClustering()
	vtkIntArray *InitialClustering;

	/// Parameter to enable Graphical Display of the clustering 
	/// 0: no display  1: display after convergence   2: display during convergence (might be slow)
	int Display;
	
	/// initialize the arrays (accumulators etc)
	virtual void Init();
	
	/// re-coputes the statistics of the clustering
	void ReComputeStatistics();	

	/// the maximum number of convergences. 
	/// When the number max of convergences is reached, the clustering is considered to be finished
	/// Note: this is only a secutity system to prevent hanging. The clustering heoretically allways finishes
	int MaxNumberOfConvergences;

	/// When the MaxNumberOfLoops is reached, the clustering is considered to have converged
	/// Note: this is only a secutity system to prevent hanging. The clustering theoretically allways converge
	int MaxNumberOfLoops;

	/// Parameter to enable Text Output
	int ConsoleOutput;
	
	/// Picks the initial samples
	void InitSamples(vtkIdList *List);

	/// a independent method to randomly initialize regions according to the weights
	virtual void ComputeInitialRandomSampling(vtkIdList *List,vtkIntArray *Sampling,int NumberOfRegions);

	/// Paramter defining th initial sampling type
	/// 1: random initialisation  2: weight-based initialisation
	int InitialSamplingType;
	
	/// the number of wanted clusters
	int NumberOfClusters;

	/// detects the clusters made of several disconnected components
	/// their number is the returned integer. for each of these clusters
	/// only the component with the biggest area is kept. The other components
	/// are reset to the NULL cluster.
	int CleanClustering();

	// assigns the Items belonging to the NULL to other clusters, with a greedy algorithm
	void FillHolesInClustering(vtkIntArray *Clustering);

	/// allocates the needed memory
	void Allocate();

	/// the queue containing the edges between two clusters
	std::queue<EdgeType> EdgeQueue;

	/// re-compute the list of edges between two different clusters (usefull after cleaning, initialization...)
	virtual void FillQueuesFromClustering();

	/// this is the method that actually minimizes the energy term by moidifying the clustering
	virtual void MinimizeEnergy();
	
	/// this methods performs one minimization loop on the boundary edges;
	virtual int ProcessOneLoop();
	
	/// the array containing the number of items inside each cluster
	vtkIntArray *ClustersSizes;

	/// array containing the last time an edge was visited (usefull to remove duplicate entries in the queue)
	unsigned char *EdgesLastLoop;

	/// value to keep a relative notion of NumberOfloops
	unsigned char RelativeNumberOfLoops;

	/// array containing the last time a cluster was modified (usefull for speed improvement)
	int *ClustersLastModification;
	
	/// clears the optimization context 
	/// (usefull when the boundary edges list was modified externally e.g. after initialization etc...)
	void SetAllClustersToModified();

	/// the number of times the edges queue has been processed (defines the "time")
	int NumberOfLoops;

	/// the frame number when Display=3
	int FrameNumber;

	/// Timer for speed measures
	vtkTimerLog *Timer;

	/// The clustering start time is stored here
	double StartTime;

	/// returns the global energy value, also computes the difference between last time it was computed
	long double ComputeGlobalEnergy();

	/// method to enable/disable the writing of a file containing the energy values during clustering
	void SetComputeAndSaveEnergy(int S) {this->ComputeAndSaveEnergy=S;};

	/// the flag defining whether or not the energy value is stored to a file
	int ComputeAndSaveEnergy;

	/// small method to write current time and energy value to the file
	void WriteToGlobalEnergyLog();

	/// the text file where the global energy evolution is written.
	// the 4 columns are : loop number , time, energy, energy difference;
	std::ofstream GlobalEnergyLog;

	// the optionnal output directory (usefull when using the batch manager)
	char* OutputDirectory;

	// this flag determines whether the clustering will be performed first with unconstrained metric and then with constrained metric (this makes the clustering faster)
	int UnconstrainedInitializationFlag;

	// a buffer List for the method AddItemRingToProcess 
	vtkIdList *EdgeList;
};


template <class Metric, class EdgeType>
vtkIntArray* vtkUniformClustering<Metric,EdgeType>::GetClustering()
{
	return (this->Clustering);
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::SetInitialClustering(vtkIntArray *Clust)
{
	this->InitialSamplingType=2;
	this->InitialClustering=Clust;
	InitialClustering->Register(this);
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::AddItemRingToProcess(vtkIdType Item)
{
	this->GetItemEdges(Item,this->EdgeList);
	for (int i=0;i<this->EdgeList->GetNumberOfIds();i++)
		this->EdgeQueue.push(this->EdgeList->GetId(i));
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::ReComputeStatistics()
{
	int	i,Value;
	int	*Size=this->ClustersSizes->GetPointer(0);
	for	(i=0;i<this->NumberOfClusters;i++)
	{
		this->MetricContext.ResetCluster(Clusters+i);
		*Size=0;
		Size++;
	}
	for	(i=0;i<this->GetNumberOfItems();i++)
	{
		Value=this->Clustering->GetValue(i);
		if ((Value>=0)&&(Value<this->NumberOfClusters))
		{
			this->MetricContext.AddItemToCluster(i,this->Clusters+Value);
			(*this->ClustersSizes->GetPointer(Value))++;
		}
	}
	for	(i=0;i<this->NumberOfClusters;i++)
	{
	
		this->MetricContext.ComputeClusterCentroid(Clusters+i);
		this->MetricContext.ComputeClusterEnergy(Clusters+i);
	}
}

template <class Metric, class EdgeType>
int	vtkUniformClustering<Metric,EdgeType>::CleanClustering()
{
	int	i,j,l;
	int	Number,Type,Size;
	vtkIdList **Clusters;
	char *Visited;
	vtkIdType *VisitedCluster;
	int	*Sizes;
	vtkIdType I1,I2;

	std::queue<int>	Queue;
	Clusters=new vtkIdList*[this->NumberOfClusters];
	Sizes=new int[this->NumberOfClusters];
	Visited=new	char[this->GetNumberOfItems()];
	VisitedCluster=new vtkIdType[this->NumberOfClusters];
	Number=0;

	vtkIdList *IList=vtkIdList::New();

	// Detect clusters that	have several connected components
	for	(i=0;i<this->GetNumberOfItems();i++)
		Visited[i]=0;
	for	(i=0;i<this->NumberOfClusters;i++)
	{
		Clusters[i]=0;
		Sizes[i]=0;
		VisitedCluster[i]=0;
	}

	for	(i=0;i<this->GetNumberOfItems();i++)
	{
		while (Queue.empty()==0)
			Queue.pop();

		if ((Visited[i]==0)&&(this->Clustering->GetValue(i)!=this->NumberOfClusters))
		{
			Size=0;
			Queue.push(i);
			Type=this->Clustering->GetValue(i);
			while (Queue.empty()==0)
			{
				I1=Queue.front();
				Queue.pop();
				if (Visited[I1]==0)
				{
					if (this->GetItemType(I1)!=0)
						Size++;
					Visited[I1]=1;

					this->GetItemNeighbours(I1,IList);

					for	(j=0;j<IList->GetNumberOfIds();j++)
					{
						I2=IList->GetId(j);
						if ((Visited[I2]==0)&&(this->Clustering->GetValue(I2)==Type))
							Queue.push(I2);
					}
				}
			}
			if (Type!=this->NumberOfClusters)
			{
				if (VisitedCluster[Type]==0)
				{
					// first connected component
					VisitedCluster[Type]=i;
					Sizes[Type]=Size;
				}
				else
				{
					// The cluster has an other	connected component
					if (Clusters[Type]==0)
					{
						Clusters[Type]=vtkIdList::New();
						Clusters[Type]->InsertNextId(VisitedCluster[Type]);
						Clusters[Type]->InsertNextId(Sizes[Type]);
					}
					Clusters[Type]->InsertNextId(i);
					Clusters[Type]->InsertNextId(Size);
				}
			}																  
		}
	}

	for	(i=0;i<this->GetNumberOfItems();i++)
		Visited[i]=0;

	// Reset the smallest unconnected components to	NULLCLUSTER
	int	Sizemax, Imax=0;
	for	(i=0;i<this->NumberOfClusters;i++)
	{
		if ((Clusters[i]!=0)&&(this->IsClusterCleanable(i)))
		{
			Number++;
			Sizemax=0;
			// Detect for each cluster,	the	biggest	connected component;
			for	(j=0;j<Clusters[i]->GetNumberOfIds()/2;j++)
			{
				if (Sizemax<Clusters[i]->GetId(2*j+1))
				{
					Sizemax=Clusters[i]->GetId(2*j+1);
					Imax=j;
				}
			}

			// Reset the smallest components to	-1
			for	(j=0;j<Clusters[i]->GetNumberOfIds()/2;j++)
			{
				if (j!=Imax)
				{
					while (Queue.empty()==0)
						Queue.pop();
					Queue.push(Clusters[i]->GetId(2*j));
					Type=this->Clustering->GetValue(Clusters[i]->GetId(2*j));
					while (Queue.empty()==0)
					{
						I1=Queue.front();
						Queue.pop();
						if (Visited[I1]==0)
						{
							Visited[I1]=1;
							this->Clustering->SetValue(I1,this->NumberOfClusters);
							this->GetItemNeighbours(I1,IList);
							for	(l=0;l<IList->GetNumberOfIds();l++)
							{
								I2=IList->GetId(l);
								if (this->Clustering->GetValue(I2)==Type)
									Queue.push(I2);
							}
						}
					}
				}
			}
			Clusters[i]->Delete();
		}
	}

	delete [] Clusters;
	delete [] Visited;
	delete [] VisitedCluster;
	delete [] Sizes;
	IList->Delete();
	
	return (Number);
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::FillHolesInClustering(vtkIntArray *Clustering)
{
	std::queue<int>	IQueue;
	vtkIdType	i;
	int	Cluster1,Cluster2;
	vtkIdType I1,I2;
	int	Edge;
	int	TempValue;
	int	InitialNumberOfProblems=0,NumberOfProblems;

	for	(i=0;i<this->GetNumberOfItems();i++)
	{
		Cluster1=Clustering->GetValue(i);
		if ((Cluster1<0)||(Cluster1>=this->NumberOfClusters))
			InitialNumberOfProblems++;
	}

	for (i=0;i<this->GetNumberOfEdges();i++)
	{
		this->GetEdgeItems(i,I1,I2);
		if ((I2>=0)&&(I1>=0))
		{
			Cluster1=Clustering->GetValue(I1);		
			Cluster2=Clustering->GetValue(I2);	
			if ((Cluster1<0)||(Cluster1>=this->NumberOfClusters))
			{
				if ((Cluster2>=0)&&(Cluster2<this->NumberOfClusters))
					IQueue.push(i);
			}
			else
			{
				if ((Cluster2<0)||(Cluster2>=this->NumberOfClusters))
					IQueue.push(i);
			}
		}
	}

	vtkIdList *EList=vtkIdList::New();

	while (IQueue.size()!=0)
	{
		Edge=IQueue.front();
		IQueue.pop();
		this->GetEdgeItems(Edge,I1,I2);
		if ((I1>=0)&&(I2>=0))
		{
			Cluster1=Clustering->GetValue(I1);
			Cluster2=Clustering->GetValue(I2);

			if ((Cluster1==this->NumberOfClusters))
			{
				TempValue=I1;
				I1=I2;
				I2=TempValue;
				Cluster1=Cluster2;
				Cluster2=this->NumberOfClusters;
			}
			
			if ((Cluster1!=this->NumberOfClusters)&&(Cluster2==this->NumberOfClusters)
					&&(this->ConnexityConstraintProblem(I2,Edge,Cluster2,Cluster1)==0))
			{
				Clustering->SetValue(I2,Cluster1);
				this->GetItemEdges(I2,EList);
				int j;
				for	(j=0;j<EList->GetNumberOfIds();j++)
				{
					IQueue.push(EList->GetId(j));
				}
			}
		}
	}

	EList->Delete();
	EList=0;
	NumberOfProblems=0;

	for	(i=0;i<this->GetNumberOfItems();i++)
	{
		Cluster1=Clustering->GetValue(i);
		if ((Cluster1<0)||(Cluster1>=this->NumberOfClusters))
			NumberOfProblems++;
	}
	
	if ((this->ConsoleOutput>0)&&(NumberOfProblems>0))
		cout<<endl<<"WARNING : The number	of uncorrectly assigned	Items was reduced from "<<
		InitialNumberOfProblems<<" to "<<NumberOfProblems<<" Problems"<<endl;
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::FillQueuesFromClustering()
{
	vtkIdType Cluster1,Cluster2;
	vtkIdType I1,I2;

	while (this->EdgeQueue.size())
		this->EdgeQueue.pop();

	for	(int i=0;i<this->GetNumberOfEdges();i++)
	{
		this->GetEdgeItems(i,I1,I2);
		if (I2>=0)
		{
			Cluster1=this->Clustering->GetValue(I1);
			Cluster2=this->Clustering->GetValue(I2);
			if (Cluster1!=Cluster2)
				this->EdgeQueue.push(i);
		}
	}
	this->EdgeQueue.push(-1);
}

template <class Metric, class EdgeType>
vtkIntArray* vtkUniformClustering<Metric,EdgeType>::ProcessClustering(vtkIdList *List)
{
	vtkTimerLog	*Timer=vtkTimerLog::New();
	if (this->NumberOfClusters==0)
	{
		cout<<"Problem!!! must set NumberOfClusters	to more	than zero!"<<endl;
		return (0);
	}

	this->Init();
	this->InitSamples(List);
	
	char FileName[1000];
	if (this->OutputDirectory)
	{
		strcpy(FileName,this->OutputDirectory);
		strcat(FileName,"energy.txt");
	}
	else
		strcpy(FileName,"energy.txt");
		
	
	if (this->ComputeAndSaveEnergy)
		this->GlobalEnergyLog.open (FileName, ofstream::out | ofstream::trunc);

	this->CreateWindows();

	if (this->ConsoleOutput)
		cout<<"Clustering......"<<endl;

	if (this->UnconstrainedInitializationFlag)
	{
		if (this->ConsoleOutput)
			cout<<"Performing unconstrained initialization"<<endl;
		this->MetricContext.SetConstrainedClustering(0);
	}
	
	this->StartTime=Timer->GetUniversalTime();
	Timer->StartTimer();
	this->MinimizeEnergy();
	Timer->StopTimer();
	
	if (this->ComputeAndSaveEnergy)
	{
		this->ReComputeStatistics();
		GlobalEnergyLog<<"Final Energy :"<<setprecision(15)<<this->ComputeGlobalEnergy()<<endl;		
		this->GlobalEnergyLog.close();
	}

	if (this->ConsoleOutput!=0)
	{
		cout<<"The clustering took :"<<Timer->GetElapsedTime()<<" seconds."<<endl;
		cout<<"Number of loops:	"<<this->NumberOfLoops<<endl;
	}

	Timer->Delete();
	this->Snapshot();
	return (this->Clustering);
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::SetAllClustersToModified()
{
	for	(int i=0;i<this->NumberOfClusters;i++)
		this->ClustersLastModification[i]=this->NumberOfLoops;
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::MinimizeEnergy()
{
	int NumberOfConvergences=0;
	int NumberOfModifications;
	int NumberOfDisconnectedClusters;

	this->FillHolesInClustering(this->Clustering);
	this->FillQueuesFromClustering();
	this->ReComputeStatistics();
	this->SetAllClustersToModified();
	vtkTimerLog *Timer=vtkTimerLog::New();
	
	while (1)
	{
		Timer->StartTimer();
		this->SwapQueues();	
		NumberOfModifications=this->ProcessOneLoop();
		Timer->StopTimer();
		
		// Display the clustering if wanted
		if (this->Display>1)
			this->Snapshot();

		// Write Energy value and computing times to file if wanted
		this->WriteToGlobalEnergyLog();
		
		if (this->ConsoleOutput>1)
		{
			cout<<(char) 13;
			if (this->ConnexityConstraint==1)
				cout<<"*";
					
			cout<<"Loop "<<this->NumberOfLoops<<", duration : "<<(int) Timer->GetElapsedTime()
				<<" s., "<<NumberOfModifications<<" Modifications            "<<std::flush;
		}
		this->NumberOfLoops++;
		if (this->RelativeNumberOfLoops==255)
		{
			//	reset the EdgesLastLoop array to cope with overflow
			for (int i=0;i<this->GetNumberOfEdges();i++)
				this->EdgesLastLoop[i]=0;
			this->RelativeNumberOfLoops=1;
		}
		else
			this->RelativeNumberOfLoops++;
		
		if ((NumberOfModifications<1)||(this->NumberOfLoops>this->MaxNumberOfLoops)
		||(((NumberOfModifications<this->GetNumberOfItems()/1000)||(this->NumberOfLoops*2>this->MaxNumberOfLoops))
				&&(NumberOfConvergences==0)))
		{
			if ((this->UnconstrainedInitializationFlag)&&(NumberOfConvergences==0))
			{
				cout<<endl<<"Unconstrained initialization done"<<endl;
				this->MetricContext.SetConstrainedClustering(1);
			}
			else
			{
				if ((NumberOfModifications)&&(this->ConsoleOutput))
					cout<<endl<<"Trigerring early convergence for speed increase";
			}
			this->ConnexityConstraint=1;		
			NumberOfConvergences++;

			NumberOfDisconnectedClusters=this->CleanClustering();
			this->FillHolesInClustering(this->Clustering);

			if (this->ConsoleOutput)
				cout<<endl<<"Convergence: "<<NumberOfDisconnectedClusters<<" disconnected classes	  "<<endl;

			if ((NumberOfDisconnectedClusters==0)&&(NumberOfModifications==0))
				break;

			if (this->NumberOfLoops>=this->MaxNumberOfLoops)
			{
				if (this->ConsoleOutput)
					cout<<"Maximum allowed number of loops reached, exiting minimization"<<endl;
				break;
			}
			
			if(NumberOfConvergences>=MaxNumberOfConvergences)
			{
				if (this->ConsoleOutput)
					cout<<"Maximum allowed number of convergences reached, exiting minimization"<<endl;
				break;
			}
			
			this->ReComputeStatistics();			
			this->FillQueuesFromClustering();
			this->SetAllClustersToModified();
		}	
	}
	Timer->Delete();
}

template <class Metric, class EdgeType>
int vtkUniformClustering<Metric,EdgeType>::ProcessOneLoop()
{
	vtkIdType Edge,I1,I2;
	int	Val1,Val2,*Size1,*Size2;

	typename Metric::Cluster *Cluster1,*Cluster2,*Cluster21,*Cluster22,*Cluster31,*Cluster32;

	Cluster21=new typename Metric::Cluster;
	Cluster22=new typename Metric::Cluster;
	Cluster31=new typename Metric::Cluster;
	Cluster32=new typename Metric::Cluster;

	// Those variables will contain the energy values for the three possible cases.
	// The "volatile" statement is here to fix some numerical issues
	// (see : http://gcc.gnu.org/bugzilla/show_bug.cgi?id=323)
	volatile double Try1,Try2,Try3;
	volatile double Try11,Try12,Try21,Try22,Try31,Try32;

	int NumberOfModifications=0;
	while (1)
	{		
		Edge=this->EdgeQueue.front();
		this->EdgeQueue.pop();
		
		if (Edge==-1)
			break;

		this->GetEdgeItems(Edge,I1,I2);

		// Check if	this edge was not already visited. 
		if ((this->EdgesLastLoop[Edge]!=this->RelativeNumberOfLoops)&&(I2>=0))
		{
			this->EdgesLastLoop[Edge]=this->RelativeNumberOfLoops;
			{
				Val1=this->Clustering->GetValue(I1);
				Val2=this->Clustering->GetValue(I2);

				if (Val1!=Val2)
				{
					if (Val1==NumberOfClusters)
					{
						// I1 is not associated. Give it to the same cluster as I2
						this->MetricContext.AddItemToCluster(I1,this->Clusters+Val2);
						this->MetricContext.ComputeClusterCentroid(this->Clusters+Val2);
						this->MetricContext.ComputeClusterEnergy(this->Clusters+Val2);
						(*this->ClustersSizes->GetPointer(Val2))++;
						this->AddItemRingToProcess(I1);
						NumberOfModifications++;
						this->Clustering->SetValue(I1,Val2);
						this->ClustersLastModification[Val2]=this->NumberOfLoops;
					}
					else
					{
						if (Val2==NumberOfClusters)
						{
							// I2 is not associated. Give it to the same cluster as I1
							this->MetricContext.AddItemToCluster(I2,this->Clusters+Val1);
							this->MetricContext.ComputeClusterCentroid(this->Clusters+Val1);
							this->MetricContext.ComputeClusterEnergy(this->Clusters+Val1);
							(*this->ClustersSizes->GetPointer(Val1))++;
							this->AddItemRingToProcess(I2);
							NumberOfModifications++;
							this->Clustering->SetValue(I2,Val1);
							this->ClustersLastModification[Val1]=this->NumberOfLoops;
						}
						else
						{
							int Result=1;

							// determine whether one of	the	two	adjacent clusters was modified.
							//	If not,	the	test is	useless, and the speed improved	:)
							if ((this->ClustersLastModification[Val1]>=this->NumberOfLoops-1)
								||(this->ClustersLastModification[Val2]>=this->NumberOfLoops-1))
							{
								Cluster1=this->Clusters+Val1;
								Cluster2=this->Clusters+Val2;

								Size1=this->ClustersSizes->GetPointer(Val1);
								Size2=this->ClustersSizes->GetPointer(Val2);

								// Compute the initial energy
								Try11=this->MetricContext.GetClusterEnergy(Cluster1);
								Try12=this->MetricContext.GetClusterEnergy(Cluster2);
								Try1=Try11+Try12;

								// Compute the energy when setting I1 to the same cluster as I2;
								if ((*Size1==1)||(this->ConnexityConstraintProblem(I1,Edge,Val1,Val2)==1))
									Try2=100000000.0;
								else
								{
									this->MetricContext.Sub(Cluster1,I1,Cluster21);
									this->MetricContext.Add(Cluster2,I1,Cluster22);
									this->MetricContext.ComputeClusterCentroid(Cluster21);
									this->MetricContext.ComputeClusterCentroid(Cluster22);
									this->MetricContext.ComputeClusterEnergy(Cluster21);
									this->MetricContext.ComputeClusterEnergy(Cluster22);
									Try21=this->MetricContext.GetClusterEnergy(Cluster21);
									Try22=this->MetricContext.GetClusterEnergy(Cluster22);
									Try2=Try21+Try22;
								}

								// Compute the energy when setting I2 to the same cluster as I1;
								if ((*Size2==1)||(this->ConnexityConstraintProblem(I2,Edge,Val2,Val1)==1))
									Try3=1000000000.0;
								else
								{
									this->MetricContext.Sub(Cluster2,I2,Cluster32);
									this->MetricContext.Add(Cluster1,I2,Cluster31);
									this->MetricContext.ComputeClusterCentroid(Cluster31);
									this->MetricContext.ComputeClusterCentroid(Cluster32);										
									this->MetricContext.ComputeClusterEnergy(Cluster31);
									this->MetricContext.ComputeClusterEnergy(Cluster32);

									Try31=this->MetricContext.GetClusterEnergy(Cluster31);
									Try32=this->MetricContext.GetClusterEnergy(Cluster32);
									Try3=Try31+Try32;

								}

								if ((Try1 <= Try2) && (Try1 <= Try3))
									Result = 1;
								else if ((Try2 < Try1) && (Try2 < Try3))
									Result = 2;
								else
									Result = 3;

								switch (Result)
								{
								case (1):
									//Don't	do anything!
									this->EdgeQueue.push(Edge);
									break;

								case (2):
									// Set I1 in the same cluster as I2
									this->Clustering->SetValue(I1,Val2);
									(*Size2)++;
									(*Size1)--;
									this->MetricContext.DeepCopy(Cluster21,Cluster1);
									this->MetricContext.DeepCopy(Cluster22,Cluster2);
									this->AddItemRingToProcess(I1);
									NumberOfModifications++;
									this->ClustersLastModification[Val1]=this->NumberOfLoops;
									this->ClustersLastModification[Val2]=this->NumberOfLoops;
									break;

								case(3):
									// Set I2 in the same cluster as I1
									this->Clustering->SetValue(I2,Val1);
									(*Size1)++;
									(*Size2)--;
									this->MetricContext.DeepCopy(Cluster31,Cluster1);
									this->MetricContext.DeepCopy(Cluster32,Cluster2);
									this->AddItemRingToProcess(I2);
									NumberOfModifications++;
									this->ClustersLastModification[Val1]=this->NumberOfLoops;
									this->ClustersLastModification[Val2]=this->NumberOfLoops;
									break;
								}
							}
							else
							{
								//Don't	do anything!
								this->EdgeQueue.push(Edge);
							}							
						}
					}
				}
			}
		}
	}
	this->EdgeQueue.push(-1);
	delete Cluster21;
	delete Cluster22;
	delete Cluster31;
	delete Cluster32;
	return (NumberOfModifications);
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::InitSamples(vtkIdList *List)
{
	vtkMath	*CRan=vtkMath::New();
	CRan->RandomSeed(5000);
	CRan->Random();
	int	i;
	int	ok;
	float test;
	int	number;

	// reset all items to the Null cluster
	for	(i=0;i<this->GetNumberOfItems();i++)
		this->Clustering->SetValue(i,NumberOfClusters);

	switch (this->InitialSamplingType)
	{
	case 0:
		// Randomly pick one item for each cluster
		for	(i=0;i<this->NumberOfClusters;i++)
		{
			ok=0;
			while (ok==0)
			{
				if (List)
					test=CRan->Random(0.0,(float) List->GetNumberOfIds()-1);
				else 
					test=CRan->Random(0.0,(float)this->GetNumberOfItems());

				number=(int)floor(test+0.5);

				if (List)
					number=List->GetId(number);			
				if (this->Clustering->GetValue(number)==NumberOfClusters)
				{
					ok=1;
					this->Clustering->SetValue(number,i);
				}
			}
		}
		break;
	case 1:
		this->ComputeInitialRandomSampling(List,Clustering,this->NumberOfClusters);
		break;
	case 2:
		for (i=0;i<this->GetNumberOfItems();i++)
			this->Clustering->SetValue(i,this->InitialClustering->GetValue(i));
	}
	CRan->Delete();
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::ComputeInitialRandomSampling(vtkIdList *List,vtkIntArray	*Sampling,int NumberOfRegions)
{
	for	(int i=0;i<this->GetNumberOfItems();i++)
		Sampling->SetValue(i,NumberOfClusters);

	int	*Items=new int[this->GetNumberOfItems()];
	int	NumberOfRemainingItems=this->GetNumberOfItems();
	int	NumberOfRemainingRegions=NumberOfRegions;
	bool Found;
	vtkIdList *IList=vtkIdList::New();
	std::queue<int>	IQueue;

	// shuffle the Ids ordering
	for	(int i=0;i<NumberOfRemainingItems;i++)
		Items[i]=i;
	std::random_shuffle(Items, Items+NumberOfRemainingItems);

	// compute total weight
	double SWeights=0;
	for	(int i=0;i<NumberOfRemainingItems;i++)
		SWeights+=this->MetricContext.GetItemWeight(i);

	// compute desired average weight
	double Weight=SWeights/((double)NumberOfRemainingRegions);

	int FirstItem=0;
	while((NumberOfRemainingItems>0)&&(NumberOfRemainingRegions>0))
	{
		Found=false;
		int Item;
		while((!Found)&&(FirstItem<this->GetNumberOfItems()))
		{
			Item=Items[FirstItem];
			if (Sampling->GetValue(Item)==NumberOfClusters)
				Found=true;
			else
				FirstItem++;
		}
		if (Found)
		{
			while(IQueue.size())
				IQueue.pop();
			IQueue.push(Item);
			SWeights=0;

			NumberOfRemainingRegions--;
			while (IQueue.size())
			{
				Item=IQueue.front();
				IQueue.pop();
				if (Sampling->GetValue(Item)==NumberOfClusters)
				{
					Sampling->SetValue(Item,NumberOfRemainingRegions);
					SWeights+=this->MetricContext.GetItemWeight(Item);
					NumberOfRemainingItems--;
					this->GetItemNeighbours(Item,IList);

					for	(int j=0;j<IList->GetNumberOfIds();j++)
						IQueue.push(IList->GetId(j));

					if (SWeights>Weight)
						break;
				}
			}
		}
	}
	delete [] Items;
	IList->Delete();

	if (NumberOfRemainingRegions==0)
		return;

	vtkMath	*CRan=vtkMath::New();
	CRan->RandomSeed(5000);
	CRan->Random();
	while(NumberOfRemainingRegions>0)
	{
		Found=false;
		while (!Found)
		{
			float test=CRan->Random(0.0,(float)this->GetNumberOfItems()-1);
			int number=(int) floor(((double)test)+0.5);
			if (this->ClustersSizes->GetValue(this->Clustering->GetValue(number))>1)
			{
				NumberOfRemainingRegions--;
				Found=true;
				Sampling->SetValue(number,NumberOfRemainingRegions);					
			}
		}
	}
	CRan->Delete();
}

template <class Metric, class EdgeType>
long double vtkUniformClustering<Metric,EdgeType>::ComputeGlobalEnergy()
{
	// Global Energy is computed in two steps to improve precision
	int i;
	long double Energy=0;
	int NumberOfBins=(int) sqrt((double)this->NumberOfClusters);
	int Bin=0;
	int BinCount=0;
	long double *Bins=new long double [NumberOfBins];
	
	for (i=0;i<NumberOfBins;i++)
		Bins[i]=0;
		
	for (i=0;i<this->NumberOfClusters;i++)
	{
		Bins[Bin]+=(long double) this->MetricContext.GetClusterEnergy(this->Clusters+i);
		BinCount++;
		if (BinCount==NumberOfBins)
		{
			if (Bin<NumberOfBins-1)
				Bin++;
		}
	}
	
	for (i=0;i<NumberOfBins;i++)
		Energy+=Bins[i];

	delete [] Bins;
	return (Energy);
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::Allocate()
{
	this->ClustersSizes=vtkIntArray::New();
	this->ClustersSizes->SetNumberOfValues(this->NumberOfClusters);

	this->Clustering=vtkIntArray::New();
	this->Clustering->SetNumberOfValues(this->GetNumberOfItems());

	this->ClustersLastModification=new int [this->NumberOfClusters];

	this->EdgesLastLoop=new	unsigned char[this->GetNumberOfEdges()];
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::Init()
{
	int	i;

	this->Allocate();
	this->BuildMetric();

	for	(i=0;i<this->NumberOfClusters;i++)
		this->ClustersSizes->SetValue(i,0);

	while (this->EdgeQueue.size())
		this->EdgeQueue.pop();

	for	(i=0;i<this->GetNumberOfEdges();i++)
		this->EdgesLastLoop[i]=0;

	this->RelativeNumberOfLoops=1;
	this->NumberOfLoops=0;
}

template <class Metric, class EdgeType>
void vtkUniformClustering<Metric,EdgeType>::WriteToGlobalEnergyLog()
{
	if (this->ComputeAndSaveEnergy)
	{
		double Time;
		Time=Timer->GetUniversalTime();
		GlobalEnergyLog<<this->NumberOfLoops<<" "<<Time-this->StartTime
		<<" "<<setprecision(15)<<this->ComputeGlobalEnergy()<<setprecision(6)<<endl;
	}
}

template <class Metric, class EdgeType>
vtkUniformClustering<Metric,EdgeType>::vtkUniformClustering()
{
	this->Clustering=0;
	this->InitialClustering=0;
	this->ClustersSizes=0;
	this->EdgesLastLoop=0;
	this->ClustersLastModification=0;
	this->Display=0;
	this->ConsoleOutput=0;
	this->MaxNumberOfConvergences=1000000000;
	this->MaxNumberOfLoops=500;
	this->NumberOfClusters=0;
	this->InitialSamplingType=1;
	this->ConnexityConstraint=0;
	this->ComputeAndSaveEnergy=0;
	this->OutputDirectory=0;
	this->Clusters=0;
	this->UnconstrainedInitializationFlag=0;
	this->EdgeList=vtkIdList::New();
}

template <class Metric, class EdgeType>
vtkUniformClustering<Metric,EdgeType>::~vtkUniformClustering()
{
	if (this->Clustering)
		this->Clustering->Delete();

	if (this->InitialClustering)
		this->InitialClustering->Delete();

	if (this->ClustersSizes)
		this->ClustersSizes->Delete();

	if (this->EdgesLastLoop)
		delete [] this->EdgesLastLoop;

	if (this->ClustersLastModification)
		delete [] this->ClustersLastModification;

	while (this->EdgeQueue.size())
		this->EdgeQueue.pop();

	if (this->Clusters)
		delete [] this->Clusters;

	this->EdgeList->Delete();
}
#endif
