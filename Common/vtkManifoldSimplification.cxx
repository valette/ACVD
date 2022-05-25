#include <vtkObjectFactory.h>
#include <vtkQuadric.h>
#include "vtkQuadricTools.h"
#include "vtkSurfaceIterators.h"
#include "vtkManifoldSimplification.h"


void vtkManifoldSimplification::Simplify()
{
	this->AllocateMemory();

	// Compute quadrics for each vertex
	vtkQuadricTools *Tool=vtkQuadricTools::New();	
	for (vtkIdType Point=0;Point!=this->Input->GetNumberOfPoints();Point++)
		Tool->GetPointQuadric(this->Input,Point,this->Quadrics[Point]);
	Tool->Delete();

	// fill up queue
	for (vtkIdType Edge=0;Edge!=this->Input->GetNumberOfEdges();Edge++)
		this->UpdateEdgePriority(Edge);

	vtkIdList *EdgesToUpdate=vtkIdList::New();
	vtkIdList *Neighbours1=vtkIdList::New();
	vtkIdList *Neighbours2=vtkIdList::New();
	vtkSurfaceVertexRingRandomIterator Iterator;
	Iterator.SetInputData(this->Input);

	int NextNumberOfPoints=this->Input->GetNumberOfPoints();
	int CurrentNumberOfPoints=this->Input->GetNumberOfPoints();
	
	// remove disconnected points from accounting
	for (vtkIdType Vertex=0;Vertex!=this->Input->GetNumberOfPoints();Vertex++)
	{
		if (this->Input->GetValence(Vertex)==0)
			CurrentNumberOfPoints--;
	}

	if (this->Display)
	{
		this->Window=RenderWindow::New();
		this->Window->SetInputData(this->Input);
		this->Window->SetWindowName("Simplification");	
	}

	// iteratively collapse the edges
	while ((this->EdgesQueue->GetNumberOfItems()!=0)
		&&(CurrentNumberOfPoints>this->NumberOfOutputVertices))
	{
		if ((this->Display!=0)&&(CurrentNumberOfPoints==NextNumberOfPoints))
		{
			this->Input->DisplayMeshProperties();
			if (this->Window->GetEdgesActor())
				Window->DisplayInputEdges();
			this->Window->Render();
			this->Window->Interact();
//			NextNumberOfPoints-=10;
			NextNumberOfPoints=(NextNumberOfPoints*7)/10;
		}

		double Priority;
		vtkIdType Edge=this->EdgesQueue->Pop(0,Priority);
		int Direction=Edge & 1;
		Edge=Edge/2;

		// first test whether the edge is contractible or not topological constraint
		vtkIdType v1,v2;
		this->Input->GetEdgeVertices(Edge,v1,v2);
//		cout<<"Collapsing Edge "<<Edge<<" with vertices "<<v1<<" and "<<v2<<" , priority="<<Priority<<endl;

		EdgesToUpdate->Reset();
		Neighbours1->Reset();
		Neighbours2->Reset();
		
		Iterator.InitTraversal(v1);
		vtkIdType Vertex=Iterator.GetNextVertex();
		while (Vertex!=-1)
		{
			Neighbours1->InsertNextId(Vertex);
			EdgesToUpdate->InsertUniqueId(Iterator.GetEdge());
			Vertex=Iterator.GetNextVertex();
		}

		Iterator.InitTraversal(v2);
		Vertex=Iterator.GetNextVertex();
		while (Vertex!=-1)
		{
			Neighbours2->InsertNextId(Vertex);
			EdgesToUpdate->InsertUniqueId(Iterator.GetEdge());
			Vertex=Iterator.GetNextVertex();
		}
		
		Neighbours1->IntersectWith(Neighbours2);
		
		if (Neighbours1->GetNumberOfIds()==this->Input->GetEdgeNumberOfAdjacentFaces(Edge))
		{
			// Check whether some triangles orientation will be flipped...
			bool OK=true;
			vtkIdType f1,f2;
			double Point1[3],Point2[3],Point3[3],Point4[3];
			vtkIdType DisappearingVertex;
			
			if (Direction==0)
			{
				DisappearingVertex=v2;
				this->Input->GetVertexNeighbourFaces(v2,FacesList);
				this->Input->GetPoint(v1,Point1);
			}
			else
			{
				DisappearingVertex=v1;
				this->Input->GetVertexNeighbourFaces(v1,FacesList);
				this->Input->GetPoint(v2,Point1);
			}

			this->Input->GetEdgeFaces(Edge,f1,f2);
			FacesList->DeleteId(f1);
			FacesList->DeleteId(f2);

			for (int i=0;i!=FacesList->GetNumberOfIds();i++)
			{
				vtkIdType Vertices[3];
				double Normal1[3]={0,0,0};
				double Normal2[3]={0,0,0};
				this->Input->GetFaceVertices(FacesList->GetId(i),Vertices[0],Vertices[1],Vertices[2]);
				this->Input->GetPoint(Vertices[0],Point2);
				this->Input->GetPoint(Vertices[1],Point3);
				this->Input->GetPoint(Vertices[2],Point4);
				vtkTriangle::ComputeNormalDirection(Point2,Point3,Point4,Normal1);
				for (int i=0;i<3;i++)
				{
					if (Vertices[i]==DisappearingVertex)
					{
						switch(i)
						{
						case 0:
							vtkTriangle::ComputeNormalDirection(Point1,Point3,Point4,Normal2);
							break;
						case 1:
							vtkTriangle::ComputeNormalDirection(Point2,Point1,Point4,Normal2);
							break;
						case 2:
						default:
							vtkTriangle::ComputeNormalDirection(Point2,Point3,Point1,Normal2);
							break;
						}
						break;
					}
				}
				if (vtkMath::Dot(Normal1,Normal2)<0)
				{
					OK=false;
					break;
				}
			}
			
			if (OK)
			{
				double *NewQuadric;
				// merge the vertices
				if (Direction==0)
				{
					this->Input->MergeVertices(v1,v2);
					NewQuadric=this->Quadrics[v1];
				}
				else
				{
					this->Input->MergeVertices(v2,v1);
					NewQuadric=this->Quadrics[v2];
				}

				// update surrounding edges in the queue
				for (int i=0;i!=EdgesToUpdate->GetNumberOfIds();i++)
					this->UpdateEdgePriority(EdgesToUpdate->GetId(i));

				CurrentNumberOfPoints--;

				for (int i=0;i!=10;i++)
					NewQuadric[i]=this->Quadrics[v1][i]+this->Quadrics[v2][i];
			}
			else
			{
				this->NonContractibleEdges->SetValue(Edge,1);
				this->UpdateEdgePriority(Edge);
			}
		}
		else
		{
			this->NonContractibleEdges->SetValue(Edge,1);
			this->UpdateEdgePriority(Edge);
		}
	}
	
	if (this->Display!=0)
	{
		this->Input->DisplayMeshProperties();
		this->Window->SetWindowName("Final Simplification");
		if (this->Window->GetEdgesActor())
			Window->DisplayInputEdges();
		this->Window->Render();
		this->Window->Interact();
	}
	
	Neighbours1->Delete();
	Neighbours2->Delete();
	EdgesToUpdate->Delete();
	this->ReleaseMemory();
}

void vtkManifoldSimplification::UpdateEdgePriority(vtkIdType Edge)
{
	this->EdgesQueue->DeleteId(Edge*2);
	this->EdgesQueue->DeleteId(Edge*2+1);

	if ((this->NonContractibleEdges->GetValue(Edge)==1)
		||(this->Input->IsEdgeActive(Edge)==0))
		{
			this->NonContractibleEdges->SetValue(Edge,0);
			return;
		}

	vtkIdType v1,v2;

	this->Input->GetEdgeVertices(Edge,v1,v2);

	double Point[3];
	this->Input->GetPoint(v1,Point);
	double CurrentError=vtkQuadricTools::Evaluate(this->Quadrics[v1],Point);
	this->Input->GetPoint(v2,Point);
	
//	cout<<"Error1 = "<<CurrentError;
	double Error2=vtkQuadricTools::Evaluate(this->Quadrics[v2],Point);
//	cout<<" . Error 2= "<<Error2;
	CurrentError+=Error2;
	
	double Quadric[10];
	for (int i=0;i!=10;i++)
		Quadric[i]=this->Quadrics[v1][i]+this->Quadrics[v2][i];

	for (int Direction=0;Direction!=2;Direction++)
	{
		if (Direction==0)
			this->Input->GetPoint(v1,Point);
		else
			this->Input->GetPoint(v2,Point);

		double NewError=vtkQuadricTools::Evaluate(Quadric,Point);
		this->EdgesQueue->Insert(NewError-CurrentError,Edge*2+Direction);
//		cout<<" Error direction "<<Direction<<"="<<NewError<<" .";
	}
//	cout<<endl;
}


void vtkManifoldSimplification::AllocateMemory()
{
	int NumberOfVertices=this->Input->GetNumberOfPoints();

	this->QuadricsArray=new double[10*NumberOfVertices];
	this->Quadrics=new double*[NumberOfVertices];

	for (int i=0;i<NumberOfVertices;i++)
		this->Quadrics[i]=this->QuadricsArray+10*i;

	this->EdgesQueue->Allocate(this->Input->GetNumberOfEdges()*2);
	this->NonContractibleEdges->SetNumberOfValues(this->Input->GetNumberOfEdges());
	
	for (vtkIdType Edge=0;Edge!=this->Input->GetNumberOfEdges();Edge++)
		this->NonContractibleEdges->SetValue(Edge,0);
}

void vtkManifoldSimplification::ReleaseMemory()
{
	delete [] this->QuadricsArray;
	delete [] this->Quadrics;
}

vtkStandardNewMacro(vtkManifoldSimplification);

vtkManifoldSimplification::vtkManifoldSimplification()
{
	this->Input=0;
	this->Display=0;
	this->NumberOfOutputVertices=100;
	this->EdgesQueue=vtkPriorityQueue::New();
	this->NonContractibleEdges=vtkBitArray::New();
	this->FacesList=vtkIdList::New();
}

vtkManifoldSimplification::~vtkManifoldSimplification()
{
	this->EdgesQueue->Delete();

	if (this->Input)
		this->Input->UnRegister(this);

	this->NonContractibleEdges->Delete();
	this->FacesList->Delete();
}
