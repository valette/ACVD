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

#include <math.h>
#include <sstream>
#include <vtkObjectFactory.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkMapper.h>
#include <vtkVRMLImporter.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPriorityQueue.h>
#include <vtkEdgeTable.h>
#include <vtkOBJReader.h>
#include <vtkSTLReader.h>
#include <vtkMeshQuality.h>

#include "vtkSurface.h"
#include "vtkOFFReader.h"
#include "vtkVolumeProperties.h"

// this variable allows the creation of render windows when cleaning a mesh.
//#define DISPLAYMESHCLEANING

#ifdef DISPLAYMESHCLEANING
#include "RenderWindow.h"
#endif

bool sortByDecreasingSecond( const std::pair< int,int > &a, const std::pair< int,int > &b ) {
    return ( a.second > b.second );
}

vtkSurface *vtkSurface::GetBiggestConnectedComponents( int numberOfComponents )
{
	vtkIdListCollection *Collection=this->GetConnectedComponents();
	std::vector< std::pair< int, int > > components;

	if ( Collection->GetNumberOfItems() < numberOfComponents )
		numberOfComponents = Collection->GetNumberOfItems();

	for (int i=0;i<Collection->GetNumberOfItems();i++)
		components.push_back(
			std::make_pair( i, Collection->GetItem(i)->GetNumberOfIds() ) );

    sort(components.begin(), components.end(), sortByDecreasingSecond); 


	for ( int i = 0; i < components.size(); i++) {
		cout<< components[ i ].first << ":" << components[ i ].second << endl;
	}

	vtkIdType *Ids=new vtkIdType[this->GetNumberOfPoints()];

	for (vtkIdType i=0;i<this->GetNumberOfPoints();i++)
		Ids[i]=-1;

	vtkSurface *NewMesh=vtkSurface::New();

	for ( int comp = 0; comp < numberOfComponents; comp++ ) {

		vtkIdList *List=Collection->GetItem( components[ comp ].first );
		double Point[ 3 ];
		for (vtkIdType i=0;i<List->GetNumberOfIds();i++) {

			vtkIdType PtId=List->GetId(i);
			this->GetPoint(PtId,Point);
			Ids[PtId] = NewMesh->AddVertex(Point);

		}

	}

	int NumCells=this->GetNumberOfCells();
	vtkIdType *pts;
	vtkIdType npts, NewVertices[1000];

	for (vtkIdType i=0;i<NumCells;i++)
	{
		this->GetFaceVertices(i,npts, pts);
		bool keep=true;
		for (int j=0;j<npts;j++)
			if (Ids[pts[j]]<0)
				keep=false;

		if (keep)
		{
			for (int j=0;j<npts;j++)
			{
				vtkIdType NewVertex=Ids[pts[j]];
				if (NewVertex<0)
				{
					cout<<"Error in GetBiggestConnectedComponent()"<<endl;
					exit(1);
				}
				NewVertices[j]=NewVertex;
			}
			NewMesh->AddPolygon(npts,NewVertices);
		}
	}

	delete [] Ids;
	return (NewMesh);
}


vtkSurface *vtkSurface::GetBiggestConnectedComponent()
{
	vtkIdListCollection *Collection=this->GetConnectedComponents();

	// get the biggest component vertices
	vtkIdType Id;
	int Max=0;
	for (int i=0;i<Collection->GetNumberOfItems();i++)
	{
		int Size=Collection->GetItem(i)->GetNumberOfIds();
		if (Size>Max)
		{
			Max=Size;
			Id=i;
		}
	}

	vtkIdList *List=Collection->GetItem(Id);

	vtkIdType *Ids=new vtkIdType[this->GetNumberOfPoints()];

	for (vtkIdType i=0;i<this->GetNumberOfPoints();i++)
		Ids[i]=-1;

	vtkSurface *NewMesh=vtkSurface::New();
	double Point[3];
	for (vtkIdType i=0;i<List->GetNumberOfIds();i++)
	{
		vtkIdType PtId=List->GetId(i);
		this->GetPoint(PtId,Point);
		Ids[PtId]=NewMesh->AddVertex(Point);
	}

	int NumCells=this->GetNumberOfCells();
	vtkIdType *pts;
	vtkIdType npts, NewVertices[1000];

	for (vtkIdType i=0;i<NumCells;i++)
	{
		this->GetFaceVertices(i,npts, pts);
		bool keep=true;
		for (int j=0;j<npts;j++)
			if (Ids[pts[j]]<0)
				keep=false;

		if (keep)
		{
			for (int j=0;j<npts;j++)
			{
				vtkIdType NewVertex=Ids[pts[j]];
				if (NewVertex<0)
				{
					cout<<"Error in GetBiggestConnectedComponent()"<<endl;
					exit(1);
				}
				NewVertices[j]=NewVertex;
			}
			NewMesh->AddPolygon(npts,NewVertices);
		}
	}

	delete [] Ids;
	return (NewMesh);
}

vtkSurface* vtkSurface::CleanMemory()
{
	vtkIdType NumberOfPoints=this->GetNumberOfPoints();
	vtkIdType *Ids=new vtkIdType[NumberOfPoints];
	vtkSurface *Output=vtkSurface::New();

	// Clean the Points
	for (vtkIdType Id=0;Id<NumberOfPoints;Id++)
	{
		if (this->GetValence(Id)!=0)
		{
			double Point[3];
			this->GetPoint(Id,Point);
			Ids[Id]=Output->AddVertex(Point);
		}
	}

	// Clean the Polygons
	vtkIdType NumberOfCells=this->GetNumberOfCells();
	for (vtkIdType Id=0;Id<NumberOfCells;Id++)
	{
		if (this->IsFaceActive(Id)==1)
		{
			vtkIdType *Vertices,NVertices;
			this->GetFaceVertices(Id, NVertices,Vertices);

			// Change vertices Ids according to the cleaning phase
			for (int i=0;i<NVertices;i++)
				Vertices[i]=Ids[Vertices[i]];

			Output->AddPolygon(NVertices,Vertices);
		}
	}

	delete [] Ids;
	return (Output);
}

void vtkSurface::EnsureOutwardsNormals()
{
	vtkVolumeProperties *Volume=vtkVolumeProperties::New();
	
	Volume->SetInputData(this);
	if (Volume->GetSignedVolume()<0)
		this->SwitchOrientation();
	Volume->Delete();
}


double vtkSurface::GetEdgeLength(vtkIdType Edge)
{
	vtkIdType v1,v2;
	double P1[3],P2[3];
	this->GetEdgeVertices(Edge,v1,v2);
	this->GetPoint(v1,P1);
	this->GetPoint(v2,P2);
	return (sqrt(vtkMath::Distance2BetweenPoints(P1,P2)));
}

double TriangleMinAngle( vtkSurface *mesh, vtkIdType v1, vtkIdType v2, vtkIdType v3 )
{
  double p0[3],p1[3],p2[3];
  double a[3],b[3],c[3];
  double a2,b2,c2,alpha,beta,gamma;
  const double normal_coeff = .3183098861837906715377675267450287;

  mesh->GetPointCoordinates( v1, p0 );
  mesh->GetPointCoordinates( v2, p1 );
  mesh->GetPointCoordinates( v3, p2 );
  a[0] = p1[0]-p0[0];
  a[1] = p1[1]-p0[1];
  a[2] = p1[2]-p0[2];
 
  b[0] = p2[0]-p1[0];
  b[1] = p2[1]-p1[1];
  b[2] = p2[2]-p1[2];
 
  c[0] = p2[0]-p0[0];
  c[1] = p2[1]-p0[1];
  c[2] = p2[2]-p0[2];
 
  a2 = vtkMath::Dot(a,a);
  b2 = vtkMath::Dot(b,b);
  c2 = vtkMath::Dot(c,c);

  alpha = acos(vtkMath::Dot(b,c) / sqrt(b2 * c2));
  beta  = acos(vtkMath::Dot(c,a) / sqrt(c2 * a2));
  gamma = acos(-vtkMath::Dot(a,b) / sqrt(a2 * b2));

  alpha = alpha < beta ? alpha : beta;
  return  (alpha < gamma ? alpha : gamma) * 180. * normal_coeff;
}

void vtkSurface::WriteMeshStatisticsFile( const char* AreaFileName, const char* QFileName,
                              const char* AngleMinFileName, const char* DegreeFileName )
{
    std::ofstream area_file( AreaFileName );
    std::ofstream Q_file( QFileName );
    std::ofstream AngleMin_file( AngleMinFileName );
    std::ofstream Degree_file( DegreeFileName );

    double coeff = (12.0/sqrt(3.0));
    vtkIdType v1, v2, v3;
    double FP1[3], FP2[3], FP3[3];
    double length = 0.;
    double perimeter = 0.;
    double area = 0.;
    double lmax = 0.;

    for( vtkIdType i = 0; i < this->GetNumberOfCells( ); i++ )
    {
        this->GetFaceVertices(i,v1,v2,v3);
        this->GetPoints()->GetPoint(v1,FP1);
        this->GetPoints()->GetPoint(v2,FP2);
        this->GetPoints()->GetPoint(v3,FP3);

        length=sqrt(vtkMath::Distance2BetweenPoints(FP1,FP2));
        lmax=length;
        perimeter=length;

        length=sqrt(vtkMath::Distance2BetweenPoints(FP2,FP3));

        if( lmax < length )
            lmax=length;
        perimeter+=length;

        length=sqrt(vtkMath::Distance2BetweenPoints(FP1,FP3));
        if( lmax < length)
            lmax=length;
        perimeter+=length;
        area=vtkTriangle::TriangleArea(FP1,FP2,FP3);

        area_file <<area <<std::endl;

        Q_file <<coeff*area/(perimeter*lmax) <<std::endl;

        AngleMin_file <<TriangleMinAngle(this, v1,v2,v3) <<std::endl;
    }

    for( vtkIdType i = 0; i < this->GetNumberOfPoints( ); i++ )
        Degree_file <<this->GetValence(i) <<std::endl;

    area_file.close( );
    Q_file.close( );
    Q_file.close( );
    Degree_file.close( );
}

void vtkSurface::GetVertexNormal(vtkIdType Vertex, double *Normal)
{
	vtkIdList *FList=vtkIdList::New();
	this->GetVertexNeighbourFaces(Vertex, FList);
	vtkIdType i;
	double TempNormal[3];
	double P1[3],P2[3],P3[3];
	double Area;
	vtkIdType v1,v2,v3;
	Normal[0]=0;
	Normal[1]=0;
	Normal[2]=0;

	for (i=0;i<FList->GetNumberOfIds();i++)
	{
		this->GetFaceVertices(FList->GetId(i),v1,v2,v3);
		this->GetPointCoordinates(v1,P1);
		this->GetPointCoordinates(v2,P2);
		this->GetPointCoordinates(v3,P3);
		Area=vtkTriangle::TriangleArea(P1,P2,P3);
		vtkTriangle::ComputeNormal(P1,P2,P3,TempNormal);
		Normal[0]-=Area*TempNormal[0];
		Normal[1]-=Area*TempNormal[1];
		Normal[2]-=Area*TempNormal[2];
	}

	vtkMath::Normalize(Normal);
	FList->Delete();

}
void vtkSurface::GetTriangleNormal(vtkIdType Triangle, double *Normal)
{
	vtkIdType v1,v2,v3;
	double P1[3],P2[3],P3[3];
	this->GetFaceVertices(Triangle,v1,v2,v3);
	this->GetPointCoordinates(v1,P1);
	this->GetPointCoordinates(v2,P2);
	this->GetPointCoordinates(v3,P3);
	vtkTriangle::ComputeNormal(P1,P2,P3,Normal);
}


double vtkSurface::GetBoundingBoxDiagonalLength()
{
	double Bounds[6];
	double P1[3],P2[3];
	int i;
	this->ComputeBounds();
	this->GetBounds(Bounds);
	for (i=0;i<3;i++)
	{
		P1[i]=Bounds[2*i];
		P2[i]=Bounds[2*i+1];
	}
	return sqrt(vtkMath::Distance2BetweenPoints(P1,P2));
}

void vtkSurface::AddNoise(double Magnitude)
{
	double Strength=2*this->GetBoundingBoxDiagonalLength()*Magnitude;
	int i;
	vtkMath *Math=vtkMath::New();
	Math->RandomSeed(10);
	Math->RandomSeed(10);
	double P[3];
	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->GetPoint(i,P);
		P[0]+=(Math->Random()-0.5)*Strength;
		P[1]+=(Math->Random()-0.5)*Strength;
		P[2]+=(Math->Random()-0.5)*Strength;
		this->SetPointCoordinates(i,P);
	}
	Math->Delete();
}
void Split3(vtkSurface *Mesh, int f, int v1, int v2, int v3, int v12, int v13)
{
	Mesh->DeleteFace(f);
	Mesh->AddFace(v1,v12,v13);
	Mesh->AddFace(v12,v2,v3);
	Mesh->AddFace(v3,v13,v12);
}
void Split2(vtkSurface *Mesh, int f,int v1, int v2, int v3, int v12)
{
	Mesh->DeleteFace(f);
	Mesh->AddFace(v1,v12,v3);
	Mesh->AddFace(v12,v2,v3);
}

void vtkSurface::SplitLongEdges(double Ratio)
{
	vtkIdType i;
	vtkIdType NumberOfVertices,*Vertices;

	// prevent the vertices deletion when spliting faces
	int BackupCleanEdges=this->GetCleanEdgesState();
	int BackupCleanVertices=this->GetCleanVerticesState();
	this->SetCleanEdges(1);
	this->SetCleanVertices(0);

	// compute the average edge length multiplied by Ratio;
	vtkDoubleArray *EdgesLength=this->GetEdgeLengths();
	double Threshold=0;
	for (i=0;i<this->GetNumberOfEdges();i++)
	{
		Threshold+=EdgesLength->GetValue(i);
	}
	Threshold=Ratio*Threshold/(double) this->GetNumberOfEdges();

	vtkIdType v1,v2,v3,v12,v13,v23;
	double P1[3],P2[3],P[3];
	this->DeleteEdgeLengths();
	int Split=1;
	int *EdgesMidPoints;

	vtkIdList *FList=vtkIdList::New();
	vtkIdType j;
	int Triangles;
	while (Split==1)
	{
		Split=0;
		// determine which edges have to be split
		EdgesMidPoints=new int[this->GetNumberOfEdges()];
		for (i=0;i<this->GetNumberOfEdges();i++)
		{
			this->GetEdgeVertices(i,v1,v2);
			this->GetEdgeFaces(i,FList);

			// first check if the edge is adjacent only to triangles (and no other kind of polygons)
			Triangles=1;
			if (FList->GetNumberOfIds()==0)
				Triangles=0;
			for (j=0;j<FList->GetNumberOfIds();j++)
			{
				this->GetFaceVertices(FList->GetId(j),NumberOfVertices,Vertices);
				if (NumberOfVertices!=3)
				{
					Triangles=0;
				}
			}

			this->GetPointCoordinates(v1,P1);
			this->GetPointCoordinates(v2,P2);
			if ((Triangles==1)&&(sqrt(vtkMath::Distance2BetweenPoints(P1,P2))>Threshold))
			{
				P[0]=0.5*(P1[0]+P2[0]);
				P[1]=0.5*(P1[1]+P2[1]);
				P[2]=0.5*(P1[2]+P2[2]);
				EdgesMidPoints[i]=this->AddVertex(P);
			}
			else
				EdgesMidPoints[i]=-1;
		}

		// split the faces accordingly
		int NumCells=this->GetNumberOfCells();
		for (i=0;i<NumCells;i++)
		{
			this->GetFaceVertices(i,NumberOfVertices,Vertices);
			if (NumberOfVertices==3)
			{
				v1=Vertices[0];
				v2=Vertices[1];
				v3=Vertices[2];

				v12=EdgesMidPoints[this->IsEdge(v1,v2)];
				v13=EdgesMidPoints[this->IsEdge(v1,v3)];
				v23=EdgesMidPoints[this->IsEdge(v2,v3)];

				if (v12<0)
				{
					if (v13<0)
					{
						if (v23<0)
						{
							// leave the face unchanged
						}
						else
						{
							Split2(this,i,v2,v3,v1,v23);
							Split=1;
						}
					}
					else
					{
						if (v23<0)
						{
							Split2(this,i,v3,v1,v2,v13);
							Split=1;
						}
						else
						{
							Split3(this,i,v3,v1,v2,v13,v23);
							Split=1;
						}
					}
				}
				else
				{
					if (v13<0)
					{
						if (v23<0)
						{
							Split2(this,i,v1,v2,v3,v12);
							Split=1;
						}
						else
						{
							Split3(this,i,v2,v3,v1,v23,v12);
							Split=1;
						}
					}
					else
					{
						if (v23<0)
						{
							Split3(this,i,v1,v2,v3,v12,v13);
							Split=1;
						}
						else
						{
							this->DeleteFace(i);
							this->AddFace(v1,v12,v13);
							this->AddFace(v12,v2,v23);
							this->AddFace(v23,v3,v13);
							this->AddFace(v12,v23,v13);
							Split=1;
						}
					}
				}
			}
		}
		delete [] EdgesMidPoints;
	}
	this->SetCleanEdges(BackupCleanEdges);
	this->SetCleanVertices(BackupCleanVertices);
	this->CheckNormals();


	this->DeleteEdgeLengths();
	this->DeleteSharpVertices();
	this->DeleteTrianglesAreas();
	this->DeleteTrianglesNormals();
	this->DeleteConnectedComponents();
	this->DeleteVerticesAreas();
	this->Polys->Modified();
	this->Modified();
	FList->Delete();

}
vtkSurface *vtkSurface::Subdivide(vtkIntArray *Parent1, vtkIntArray *Parent2)
{
	vtkIntArray *EdgesChildren=vtkIntArray::New();
	vtkIdType NumEdges=this->GetNumberOfEdges();
	vtkIdType NumPoints=this->GetNumberOfPoints();
	vtkIdType NumFaces=this->GetNumberOfCells();
	vtkIdType i;

	EdgesChildren->SetNumberOfValues(NumEdges);

	vtkSurface* Output=vtkSurface::New();
	Output->Init(NumPoints+NumEdges,4*NumFaces,2*NumEdges+3*NumFaces);

	if (Parent1)
	{
		Parent1->Resize(NumPoints+NumEdges);
		Parent2->Resize(NumPoints+NumEdges);
	}

	double P[3],P1[3],P2[3];
	vtkIdType V1,V2,V3;

	// copy the old points
	for (i=0;i<NumPoints;i++)
	{
		this->GetPoint(i,P);
		Output->AddVertex(P);
	}

	// create the new points (middle of each edge)
	for (i=0;i<NumEdges;i++)
	{
		if (this->IsEdgeActive(i))
		{
			this->GetEdgeVertices(i,V1,V2);
			this->GetPoint(V1,P1);
			this->GetPoint(V2,P2);

			P[0]=0.5*(P1[0]+P2[0]);
			P[1]=0.5*(P1[1]+P2[1]);
			P[2]=0.5*(P1[2]+P2[2]);

			V3=Output->AddVertex(P);

			EdgesChildren->SetValue(i,V3);
			if (Parent1)
			{
				Parent1->SetValue(V3,V1);
				Parent2->SetValue(V3,V2);
			}
		}
	}

	// create the new cells
	int V4,V5,V6;
	for (i=0;i<NumFaces;i++)
	{
		if (this->IsFaceActive(i))
		{
			this->GetFaceVertices(i,V1,V2,V3);
			V4=EdgesChildren->GetValue(this->IsEdge(V1,V2));
			V5=EdgesChildren->GetValue(this->IsEdge(V2,V3));
			V6=EdgesChildren->GetValue(this->IsEdge(V3,V1));

			Output->AddFace(V1,V4,V6);
			Output->AddFace(V4,V2,V5);
			Output->AddFace(V5,V3,V6);
			Output->AddFace(V4,V5,V6);
		}
	}
	EdgesChildren->Delete();
	return Output;
}

void vtkSurface::SubdivideInPlace(vtkIntArray *Parent1, vtkIntArray *Parent2)
{
	vtkIntArray *EdgesChildren=vtkIntArray::New();
	vtkIdType NumEdges=this->GetNumberOfEdges();
	vtkIdType NumPoints=this->GetNumberOfPoints();
	vtkIdType NumFaces=this->GetNumberOfCells();
	vtkIdType i;

	// Set the cleaning flags
	int BackupCleanVertices=this->GetCleanVerticesState();
	int BackupCleanEdges=this->GetCleanEdgesState();
	this->SetCleanVertices(0);
	this->SetCleanEdges(1);

	EdgesChildren->SetNumberOfValues(NumEdges);

	double P[3],P1[3],P2[3];
	vtkIdType V1,V2,V3;

	if (Parent1)
	{
		Parent1->Resize(NumPoints+NumEdges);
		Parent2->Resize(NumPoints+NumEdges);
	}

	// create the new points (middle of each edge)
	for (i=0;i<NumEdges;i++)
	{
		if (this->IsEdgeActive(i))
		{
			this->GetEdgeVertices(i,V1,V2);
			this->GetPoint(V1,P1);
			this->GetPoint(V2,P2);

			P[0]=0.5*(P1[0]+P2[0]);
			P[1]=0.5*(P1[1]+P2[1]);
			P[2]=0.5*(P1[2]+P2[2]);

			V3=this->AddVertex(P);
			EdgesChildren->SetValue(i,V3);
			if (Parent1)
			{
				Parent1->SetValue(V3,V1);
				Parent2->SetValue(V3,V2);
			}
		}
	}

	// create the new cells
	int V4,V5,V6;
	for (i=0;i<NumFaces;i++)
	{
		if (this->IsFaceActive(i))
		{
			this->GetFaceVertices(i,V1,V2,V3);
			V4=EdgesChildren->GetValue(this->IsEdge(V1,V2));
			V5=EdgesChildren->GetValue(this->IsEdge(V2,V3));
			V6=EdgesChildren->GetValue(this->IsEdge(V3,V1));

			this->DeleteFace(i);
			this->AddFace(V1,V4,V6);
			this->AddFace(V4,V2,V5);
			this->AddFace(V5,V3,V6);
			this->AddFace(V4,V5,V6);
		}
	}
	EdgesChildren->Delete();
	this->Modified();

	// Restore cleaning flags
	this->SetCleanVertices(BackupCleanVertices);
	this->SetCleanEdges(BackupCleanEdges);

}


void vtkSurface::DisplayMeshProperties()
{
	std::stringstream Str;
	this->GetMeshProperties(Str);
	cout<<Str.str();
}

void vtkSurface::GetMeshProperties(std::stringstream &stream)
{
	vtkIdType i;
	stream<<"*****************************************************************************"<<endl;
	stream<<"Mesh with "
		<<this->GetNumberOfCells()<<" polygons, "
		<<this->GetNumberOfPoints()<<" points, "
		<<this->GetNumberOfEdges()<<" edges"<<endl;
	double bounds[6];
	this->GetBounds( bounds );
	stream<<"Bounding Box: [" <<bounds[0] <<", " <<bounds[2] <<", " <<bounds[4] <<"]"<<endl;
	stream<<"              [" <<bounds[1] <<", " <<bounds[3] <<", " <<bounds[5] <<"] " <<endl;
	
	vtkIdType nPoly=0,nTri=0,nQuad=0;
	vtkIdType *Vertices;
	vtkIdType NV, NumberOfEmptySlots=0;
	for (i=0;i<this->GetNumberOfCells();i++)
	{
		if (this->IsFaceActive(i)==1)
		{
			this->GetFaceVertices(i,NV,Vertices);
			if (NV==3)
				nTri++;
			else
			{
				if (NV==4)
					nQuad++;
				else
					nPoly++;
			}
		}
		else
			NumberOfEmptySlots++;

	}
	if (nPoly||nQuad)
	{
		stream<<"The mesh is made of "<<nTri<<" triangles, "<<nQuad<<" quads and "<<nPoly<<" other polygons"<<endl;
	}
	else
	{
		stream<<"The mesh is made only of "<<nTri<<" triangles"<<endl;
	}

	if (NumberOfEmptySlots>0)
		stream<<NumberOfEmptySlots<<" cells are not used (not active)"<<endl;

	int NonManifold=0;
	int Boundary=0;
	vtkIdType f1,f2;
	NumberOfEmptySlots=0;
	for (i=0;i<this->GetNumberOfEdges();i++)
	{
		if (this->IsEdgeActive(i)==1)
		{
			if (this->IsEdgeManifold(i)==0)
				NonManifold++;
			this->GetEdgeFaces(i,f1,f2);
			if (f2<0)
				Boundary++;
		}
		else
			NumberOfEmptySlots++;
	}
	stream<<NonManifold<<" non-manifold edges and "
		<<Boundary<<" boundary edges"<<endl;
	if (NumberOfEmptySlots)
		stream<<NumberOfEmptySlots<<" edges cells are not used (not active)"<<endl;

	vtkIdListCollection *Components=this->GetConnectedComponents();
	stream<<"The mesh has "<<Components->GetNumberOfItems()<<" connected components"<<endl;
    this->DeleteConnectedComponents();


	stream<<" Valences entropy: "<<this->GetValenceEntropy()<<endl;
	double a,b,c,number;
	a=0;
	b=0;
	c=0;
	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		number=this->GetValence(i);

		if (number==0)
		{
			b++;
		}
		else
		{
			a++;
			if (number!=6)
				c++;
		}
	}
	stream<<b<<" disconnected vertices, "
		<<a<<" connected vertices"<<endl;
	stream<<100.0*c/a<<" percent of irregular vertices"<<endl;

	if ((nPoly==0)&&(nQuad==0))
	{
		double Amin,Aav,Qmin,Qav,P30;
		this->ComputeTrianglesStatistics(Amin,Aav,Qmin,Qav,P30);
		stream<<"Mesh geometry quality:" <<
			endl<<"  AngleMin="<<Amin<<endl
			<<"  AverageMinAngle="<<Aav<<endl
			<<"  Qmin="<<Qmin<<endl
			<<"  Qav="<<Qav<<endl
			<<"  P30="<<P30<<endl;
	}

	stream<<"*****************************************************************************"<<endl;
}
void vtkSurface::PrintVerticesCoordinates()
{
	vtkIdType i;
	double P[3];
	cout<<"VerticesCoordinates:"<<endl;
	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->GetPointCoordinates(i,P);
		cout<<P[0]<<" "<<P[1]<<" "<<P[2]<<endl;
	}
}
void vtkSurface::PrintConnectivity()
{
	vtkIdType i;
	vtkIdType v1,v2,v3;
	cout<<"Mesh faces (v1 v2 v3)"<<endl;
	for (i=0;i<this->GetNumberOfCells();i++)
	{
		this->GetFaceVertices(i,v1,v2,v3);
		cout<<v1<<" "<<v2<<" "<<v3<<endl;
	}
}
void vtkSurface::SaveConnectivity(const char * FileName)
{
	vtkIdType  i;
	vtkIdType  v1,v2,v3;
	std::ofstream OutputFile;
	OutputFile.open (FileName, std::ofstream::out | std::ofstream::trunc);
	for (i=0;i<this->GetNumberOfCells();i++)
	{
		this->GetFaceVertices(i,v1,v2,v3);
		OutputFile<<v1<<" "<<v2<<" "<<v3<<endl;
	}
	OutputFile.close();
}

void vtkSurface::RescaleCoordinates()
{
	double Max[3],Min[3],P[3],MWidth;

	vtkIdType  i,j;
	for (i=0;i<3;i++)
	{
		Max[i]=-pow(10.0,300);
		Min[i]=pow(10.0,300);
	}
	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->GetPointCoordinates(i,P);
		for (j=0;j<3;j++)
		{
			if (Max[j]<P[j]) Max[j]=P[j];
			if (Min[j]>P[j]) Min[j]=P[j];
		}
	}

	MWidth=0;
	for (i=0;i<3;i++)
	{
		if (MWidth<Max[i]-Min[i]) MWidth=Max[i]-Min[i];
	}
	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->GetPointCoordinates(i,P);
		for (j=0;j<3;j++)
		{
			P[j]=(P[j]-Min[j])/MWidth-0.5;
		}
		this->SetPointCoordinates(i,P);
	}
}


double vtkSurface::GetDistanceBetweenVertices(vtkIdType V1,vtkIdType V2)
{
	double P1[3],P2[3];
	this->GetPointCoordinates(V1,P1);
	this->GetPointCoordinates(V2,P2);
	return (sqrt(vtkMath::Distance2BetweenPoints(P1,P2)));

}


vtkSurface* vtkSurface::CleanConnectivity(double tolerance)
{

	vtkIdList *VList1,*VList2,*VList;
	vtkIdType V1,V2,Type1,Type2,Edge;
	int NumberOfComponents;
	int NumberOfTetras;
	int i,j,k;
	vtkSurface *CleanedMesh;
	vtkPolyData *Draft=vtkPolyData::New();
	Draft->DeepCopy(this);
	vtkIdList *EdgesVertices=vtkIdList::New();
	VList=vtkIdList::New();

	double Bounds[6];
	double P1[3],P2[3];
	this->GetBounds(Bounds);
	for (i=0;i<3;i++)
	{
		P1[i]=Bounds[2*i];
		P2[i]=Bounds[2*i+1];
	}
	double BBoxDiag=sqrt(vtkMath::Distance2BetweenPoints(P1,P2));



	// Here we merge the points within the relative tolerance together
	vtkCleanPolyData *Cleaner=vtkCleanPolyData::New();
	Cleaner->SetInputData(this);
	Cleaner->SetTolerance(tolerance*BBoxDiag);
	Cleaner->Update();

	CleanedMesh=vtkSurface::New();
	CleanedMesh->CreateFromPolyData(Cleaner->GetOutput());
	Cleaner->Delete();


	// Here is the delaunay triangulation of the resulting points

	vtkDelaunay3D *Delaunay=vtkDelaunay3D::New();
	Delaunay->SetInputData(CleanedMesh);

	vtkIdListCollection *Components;
	Components=CleanedMesh->GetConnectedComponents();
	NumberOfComponents=Components->GetNumberOfItems();
	Delaunay->SetDebug(0);

	Delaunay->Update();


	vtkUnstructuredGrid *Tetras=Delaunay->GetOutput();
	NumberOfTetras=Tetras->GetNumberOfCells();
	cout<<NumberOfTetras<<" Tetras"<<endl;
	cout<<Tetras->GetNumberOfPoints()<<" Points in the tetras"<<endl;


	// Create an array containing for each vertex the index of its component;
	vtkIntArray *VerticesComponent=vtkIntArray::New();
	VerticesComponent->SetNumberOfValues(CleanedMesh->GetNumberOfPoints());
	for (i=0;i<NumberOfComponents;i++)
	{
		VList1=Components->GetItem(i);
		for (j=0;j<VList1->GetNumberOfIds();j++)
		{
			VerticesComponent->SetValue(VList1->GetId(j),i);
			V1=VList1->GetId(j);
		}
	}


	// Create a table for all the edges

	vtkEdgeTable *Edges=vtkEdgeTable::New();
	Edges->InitEdgeInsertion(CleanedMesh->GetNumberOfPoints(),1);

	for (i=0;i<Tetras->GetNumberOfCells();i++)
	{
		Tetras->GetCellPoints(i,VList);
		for (j=0;j<VList->GetNumberOfIds()-1;j++)
		{
			V1=VList->GetId(j);
			for (k=j+1;k<VList->GetNumberOfIds();k++)
			{
				V2=VList->GetId(k);
				Edge=Edges->IsEdge(V1,V2);
				if (Edge<0)
				{
					Edges->InsertEdge(V1,V2);
				}
			}
		}
	}
	cout<<Edges->GetNumberOfEdges()<<" edges for the tetraedral mesh"<<endl;


	// Push the edges into a Priority queue and store the edges vertices in EdgesVertices

	vtkPriorityQueue *Queue=vtkPriorityQueue::New();
	Queue->Allocate(Edges->GetNumberOfEdges());

	EdgesVertices->SetNumberOfIds(Edges->GetNumberOfEdges()*2);
	Edges->InitTraversal();
	Edge=Edges->GetNextEdge(V1,V2);


#ifdef 	DISPLAYMESHCLEANING
	vtkPolyData *EdgesP=vtkPolyData::New();
	EdgesP->SetPoints(this->GetPoints());
	EdgesP->Allocate(Edges->GetNumberOfEdges());
	vtkPolyData *EdgesP2=vtkPolyData::New();
	EdgesP2->SetPoints(this->GetPoints());
	EdgesP2->Allocate(Edges->GetNumberOfEdges());

	vtkIdList *vl=vtkIdList::New();
	vl->SetNumberOfIds(2);
#endif
	while (Edge>=0)
	{
		Queue->Insert(CleanedMesh->GetDistanceBetweenVertices(V1,V2),Edge);
		EdgesVertices->SetId(2*Edge,V1);
		EdgesVertices->SetId(2*Edge+1,V2);
#ifdef 	DISPLAYMESHCLEANING
		vl->SetId(0,V1);
		vl->SetId(1,V2);
		EdgesP->InsertNextCell(VTK_LINE,vl);
#endif
		Edge=Edges->GetNextEdge(V1,V2);
	}


	// pop the edges to connect unconnected components;
	while (NumberOfComponents>1)
	{
		Edge=Queue->Pop();
		if (Edge<0)
		{
			cout<<"***************PROBLEM*******************,"<<NumberOfComponents<<" remaining components"<<endl;
			break;
		}
		V1=EdgesVertices->GetId(2*Edge);
		V2=EdgesVertices->GetId(2*Edge+1);

#ifdef 	DISPLAYMESHCLEANING

		vl->SetId(0,V1);
		vl->SetId(1,V2);
		EdgesP2->InsertNextCell(VTK_LINE,vl);

#endif



		Type1=VerticesComponent->GetValue(V1);
		Type2=VerticesComponent->GetValue(V2);
		if ((CleanedMesh->IsEdge(V1,V2)<0)&&(Type1!=Type2)&&(V1!=V2))
		{
			CleanedMesh->AddEdge(V1,V2);
			VList1=Components->GetItem(Type1);
			VList2=Components->GetItem(Type2);
			for (i=0;i<VList2->GetNumberOfIds();i++)
			{
				VerticesComponent->SetValue(VList2->GetId(i),Type1);
				VList1->InsertNextId(VList2->GetId(i));
			}
			NumberOfComponents--;
		}
	}

#ifdef 	DISPLAYMESHCLEANING
	RenderWindow *Window=RenderWindow::New();
	Window->SetInputEdges(EdgesP);
	Window->SetInputData(this);
	Window->Render();
	Window->SetWindowName("1");
	Window->Interact();

	RenderWindow *Window2=RenderWindow::New();
	Window2->SetInputData(EdgesP2);
	Window2->SetInputData(this);
	Window2->Render();
	Window2->SetWindowName("2");
	Window2->AttachToRenderWindow(Window);
	Window2->Interact();


	RenderWindow *Window3=RenderWindow::New();
	Window3->SetInputData(CleanedMesh);
	Window3->DisplayInputEdges();
	Window3->Render();
	Window3->SetWindowName("3");
	Window3->AttachToRenderWindow(Window);
	Window3->Interact();

#endif

	CleanedMesh->DeleteEdgeLengths();
	CleanedMesh->DeleteConnectedComponents();

	Queue->Delete();

	Edges->Delete();
	Delaunay->Delete();
	VerticesComponent->Delete();
	cout<<"pro"<<endl;

	return (CleanedMesh);
}




// Compute the connected components of the mesh (It is based on vertices and edges
// So it works also with non manifold meshes)
vtkIdListCollection* vtkSurface::GetConnectedComponents()
{


	if (this->ConnectedComponents)
		return (this->ConnectedComponents);
	int i,j;
	vtkIdType v1,v2;


	std::queue<int> VQueue;
	int *Visited=new int [this->GetNumberOfPoints()];
	vtkIdList *Component;
	vtkIdList *VList=vtkIdList::New();
	this->ConnectedComponents=vtkIdListCollection::New();

	for (i=0;i<this->GetNumberOfPoints();i++)
		Visited[i]=0;

	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		if (Visited[i]==0)
		{
			// a new component is detected
			VQueue.push(i);
			Component=vtkIdList::New();
//			Component->Allocate(this->GetNumberOfPoints());

			while(VQueue.size())
			{
				v1=VQueue.front();
				if (Visited[v1]==0)
				{
					Component->InsertNextId(v1);

					Visited[v1]++;
					this->GetVertexNeighbours(v1,VList);
					for (j=0;j<VList->GetNumberOfIds();j++)
					{
						v2=VList->GetId(j);
						if (Visited[v2]==0)
							VQueue.push(v2);
					}
				}
				VQueue.pop();
			}


			// Compute barycenter of region, create a new node in the PolyData

			Component->Squeeze();
			this->ConnectedComponents->AddItem(Component);
			Component->Delete();
		}
	}

	VList->Delete();
	delete [] Visited;
	return (this->ConnectedComponents);
}

void vtkSurface::DeleteConnectedComponents()
{
	if (this->ConnectedComponents)
	{
		this->ConnectedComponents->Delete();
		this->ConnectedComponents=0;
	}
};



//Compute the Edges Lenghts
vtkDoubleArray* vtkSurface::GetEdgeLengths()
{
	vtkIdType i,v1,v2;
	if (this->EdgeLengths)
		return (this->EdgeLengths);

	this->EdgeLengths=vtkDoubleArray::New();
	this->EdgeLengths->SetNumberOfValues(this->GetNumberOfEdges());
	for (i=0;i<this->GetNumberOfEdges();i++)
	{
		this->GetEdgeVertices(i,v1,v2);
		this->EdgeLengths->SetValue(i,this->GetDistanceBetweenVertices(v1,v2));
	}
	return (this->EdgeLengths);
}
vtkDoubleArray* vtkSurface::GetTrianglesAreas()
{
	int i;
	vtkIdType v1,v2,v3;
	double P1[3],P2[3],P3[3];


	if (this->TrianglesAreas)
		return (this->TrianglesAreas);

	this->TrianglesAreas=vtkDoubleArray::New();
	this->TrianglesAreas->SetNumberOfValues(this->GetNumberOfCells());


	for (i=0;i<this->GetNumberOfCells();i++)
	{

		this->GetFaceVertices(i,v1,v2,v3);
		this->GetPointCoordinates(v1,P1);
		this->GetPointCoordinates(v2,P2);
		this->GetPointCoordinates(v3,P3);
		double a,b,c;
		a = vtkMath::Distance2BetweenPoints(P1,P2);
		b = vtkMath::Distance2BetweenPoints(P2,P3);
		c = vtkMath::Distance2BetweenPoints(P3,P1);
		this->TrianglesAreas->SetValue(i,(0.25* sqrt(fabs((double)4.0*a*c - (a-b+c)*(a-b+c)))));
	}

	return (this->TrianglesAreas);
}

vtkDoubleArray* vtkSurface::GetTrianglesNormals()
{
	int i;
	vtkIdType v1,v2,v3;
	double P1[3],P2[3],P3[3],N[3];


	if (this->TrianglesNormals)
		return (this->TrianglesNormals);

	this->TrianglesNormals=vtkDoubleArray::New();
	this->TrianglesNormals->SetNumberOfComponents(3);
	this->TrianglesNormals->SetNumberOfTuples(this->GetNumberOfCells());


	for (i=0;i<this->GetNumberOfCells();i++)
	{

		this->GetFaceVertices(i,v1,v2,v3);
		this->GetPointCoordinates(v1,P1);
		this->GetPointCoordinates(v2,P2);
		this->GetPointCoordinates(v3,P3);
		vtkTriangle::ComputeNormal(P1,P2,P3,N);
		this->TrianglesNormals->SetTuple(i,N);

	}

	return (this->TrianglesNormals);
}

double vtkSurface::GetFaceArea(vtkIdType Face)
{
	int j;
	vtkIdType v1,v2,v3;
	double Area;
	double Pf1[3],Pf2[3],Pf3[3];
	Area=0;
	vtkIdType *vertices, NumberOfVertices;
	this->GetFaceVertices( Face, NumberOfVertices, vertices );
	v1 = vertices[ 0 ];
	this->GetPoint(v1,Pf1);

	for (j=0;j<NumberOfVertices-2;j++)
	{
		v2=vertices[ j+1 ];
		v3=vertices[ j+2 ];
		this->GetPoint(v2,Pf2);
		this->GetPoint(v3,Pf3);
		Area+=vtkTriangle::TriangleArea(Pf1,Pf2,Pf3);
	}
	return (Area);
}

double vtkSurface::GetVertexArea(vtkIdType Vertex)
{
	vtkIdList *CList=vtkIdList::New();
	this->GetVertexNeighbourFaces(Vertex,CList);
	vtkIdType i;
	double Area;
	vtkIdType NumberOfVertices,*Vertices;

	Area=0;
	for (i=0;i<CList->GetNumberOfIds();i++)
	{
		this->GetFaceVertices(CList->GetId(i),NumberOfVertices,Vertices);
		Area+=this->GetFaceArea(CList->GetId(i))/(double) NumberOfVertices;

	}
	CList->Delete();
	return (Area);
}

vtkDoubleArray* vtkSurface::GetVerticesAreas()
{

	if (this->VerticesAreas)
		return (this->VerticesAreas);

	this->VerticesAreas=vtkDoubleArray::New();
	this->VerticesAreas->SetNumberOfValues(this->GetNumberOfPoints());

	vtkIdType i;

	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->VerticesAreas->SetValue(i,this->GetVertexArea(i));
	}

	return (this->VerticesAreas);
}

void vtkSurface::GetCellMassProperties(vtkIdType CellId, double &Area, double *Baricenter)
{
	vtkIdType  i,j;
	vtkIdType *Pts;
	vtkIdType v1,v2,v3;
	double Pf1[3],Pf2[3],Pf3[3];
	vtkIdType  NumberOfVertices;
	double Area2;

	this->GetFaceVertices(CellId,NumberOfVertices,Pts);

	Area=0;
	v1=Pts[0];
	this->GetPoint(v1,Pf1);

	Baricenter[0]=0;
	Baricenter[1]=0;
	Baricenter[2]=0;


	for (j=0;j<NumberOfVertices-2;j++)
	{
		v2=Pts[j+1];
		v3=Pts[j+2];
		this->GetPoint(v2,Pf2);
		this->GetPoint(v3,Pf3);
		Area2=vtkTriangle::TriangleArea(Pf1,Pf2,Pf3);
		for (i=0;i<3;i++)
		{
			Baricenter[i]+=Area2*(Pf1[i]+Pf2[i]+Pf3[i])/3.0;
		}
		Area+=Area2;
	}
	if (Area>0)
	{
		for (i=0;i<3;i++)
		{
			Baricenter[i]/=Area;
		}
	}
}



void vtkSurface::ComputeTrianglesStatistics(double &Amin, double & Aav,double &Qmin, double &Qav, double &P30)
{
	long double Sum,Sum2;
	double MinimumAngle,perimeter,area,Q,length,lmax;
	double FP1[3],FP2[3],FP3[3];

	int i;
	vtkIdType v1,v2,v3;

	Sum=0.0;
	Sum2=0.0;

	Amin=180;
	Qmin=10;
	double n=0.0;
	int NumberOfCells=0;

	for (i=0;i<this->GetNumberOfCells();i++)
	{
		if (this->IsFaceActive(i))
		{
			this->GetFaceVertices(i,v1,v2,v3);
			this->GetPoints()->GetPoint(v1,FP1);
			this->GetPoints()->GetPoint(v2,FP2);
			this->GetPoints()->GetPoint(v3,FP3);
			MinimumAngle=TriangleMinAngle(this, v1, v2, v3);

			if (!((MinimumAngle>=-180)&&(MinimumAngle<=180)))
				MinimumAngle=0;

			if  (MinimumAngle<Amin)
				Amin=MinimumAngle;
			Sum+=MinimumAngle;

			if (MinimumAngle<30)
				n++;

			length=sqrt(vtkMath::Distance2BetweenPoints(FP1,FP2));
			lmax=length;
			perimeter=length;

			length=sqrt(vtkMath::Distance2BetweenPoints(FP2,FP3));

			if (lmax<length) lmax=length;
			perimeter+=length;

			length=sqrt(vtkMath::Distance2BetweenPoints(FP1,FP3));
			if (lmax<length) lmax=length;
			perimeter+=length;

			area=vtkTriangle::TriangleArea(FP1,FP2,FP3);

			if (perimeter*lmax!=0)
				Q=(12.0/sqrt(3.0))*area/(perimeter*lmax);
			else
				Q=0;

			if (Qmin>Q) Qmin=Q;
			Sum2+=Q;
			NumberOfCells++;
		}
	}

	Aav=Sum/(long double) (NumberOfCells);
	Qav=Sum2/(long double) (NumberOfCells);
	P30=100.0*n/(3.0*(double) (NumberOfCells));
}


void vtkSurface::ComputeQualityHistogram(const char *FileName)
{
	double perimeter,area,Q,length,lmax;
	double FP1[3],FP2[3],FP3[3];

	int i;
	vtkIdType v1,v2,v3;

	int Histogram[100];
	for (i=0;i<100;i++)
		Histogram[i]=0;

	for (i=0;i<this->GetNumberOfCells();i++)
	{
		this->GetFaceVertices(i,v1,v2,v3);
		this->GetPoints()->GetPoint(v1,FP1);
		this->GetPoints()->GetPoint(v2,FP2);
		this->GetPoints()->GetPoint(v3,FP3);

		length=sqrt(vtkMath::Distance2BetweenPoints(FP1,FP2));
		lmax=length;
		perimeter=length;

		length=sqrt(vtkMath::Distance2BetweenPoints(FP2,FP3));

		if (lmax<length) lmax=length;
		perimeter+=length;

		length=sqrt(vtkMath::Distance2BetweenPoints(FP1,FP3));
		if (lmax<length) lmax=length;
		perimeter+=length;
		area=vtkTriangle::TriangleArea(FP1,FP2,FP3);

		Q=(12.0/sqrt(3.0))*area/(perimeter*lmax);

		Histogram[(int) floor(Q*100)]++;
	}

	std::ofstream OutputFile;
	OutputFile.open (FileName, std::ofstream::out | std::ofstream::trunc);
	for (i=0;i<100;i++)
	{
		OutputFile<<Histogram[i]<<endl;
	}
	OutputFile.close();

}
vtkPolyData * vtkSurface::GetEdgesPolyData()
{
	int i;
	vtkIdType v1,v2;
	vtkIdList *vl=vtkIdList::New();
	vl->SetNumberOfIds(2);

	vtkPolyData *Edges=vtkPolyData::New();

	Edges->SetPoints(this->GetPoints());
	Edges->Allocate(this->GetNumberOfEdges());
	for (i=0;i<this->GetNumberOfEdges();i++)
	{
		this->GetEdgeVertices(i,v1,v2);
		vl->SetId(0,v1);
		vl->SetId(1,v2);
		Edges->InsertNextCell(VTK_LINE,vl);
	}

	Edges->Modified();
	vl->Delete();
	return (Edges);
}

vtkPolyData * vtkSurface::GetVerticesPolyData()
{
	int i;
	vtkIdList *vl=vtkIdList::New();
	vl->SetNumberOfIds(1);

	vtkPolyData *Vertices=vtkPolyData::New();
	Vertices->SetPoints(this->GetPoints());
	Vertices->Allocate(this->GetNumberOfPoints());
	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		vl->SetId(0,i);
		Vertices->InsertNextCell(VTK_VERTEX,vl);
	}

	Vertices->Modified();
	vl->Delete();
	return (Vertices);
}


void vtkSurface::WriteInventor(const char *filename)
{

	std::ofstream File;
	File.open (filename, std::ofstream::out | std::ofstream::trunc);
	vtkIdType i, v1, v2, v3;
	double P[3];

	vtkPoints *points = this->GetPoints();

	File<<"#Inventor V2.1 ascii"<<endl;
	File<<"Separator {"<<endl;

	// * write out the coordinates *

	File<<"Coordinate3 {"<<endl;
	File<<"  point ["<<endl;

	for (i = 0; i < this->GetNumberOfPoints(); i++)
	{
		points->GetPoint(i,P);
		File<<P[0]<<" "<<P[1]<<" "<<P[2]<<endl;
	}

	File<<"]"<<endl;
	File<<"}"<<endl<<endl;

	// * write the faces *

	File<<"IndexedFaceSet {"<<endl;
	File<<"  coordIndex ["<<endl;

	for (i = 0; i <this->GetNumberOfCells();i++)
	{
		File<<"   ";
		this->GetFaceVertices(i,v1,v2,v3);
		File<<v1<<" "<<v2<<" "<<v3<<" -1"<<endl;
	}
	File<<"]"<<endl;
	File<<"}"<<endl<<endl;

	// * end separator *
	File<<"}"<<endl<<endl;

	File.close();
}
void vtkSurface::WriteSMF(const char *filename)
{
	std::ofstream File;
	File.open (filename, std::ofstream::out | std::ofstream::trunc);

	vtkIdType i;
	vtkIdType nverts=this->GetNumberOfPoints();
	vtkIdType nfaces=this->GetNumberOfCells();

	File << "# Generated from PLY data by ply2smf" << endl;
	File<< "# " <<  nverts<< " vertices" << endl;
	File << "# " << nfaces << " faces" << endl;

	double P[3];
	for(i=0; i<nverts; i++)
	{
		this->GetPoint(i,P);
		File << "v "<< P[0] << " " << P[1] << " " << P[2] << endl;
	}

	vtkIdType v1,v2,v3;
	for(i=0; i<nfaces; i++)
	{
		this->GetFaceVertices(i,v1,v2,v3);

		File << "f "<< v1+1 << " "<< v2+1 << " "<< v3+1 << endl;
	}
	File.close();
}


// ****************************************************************
// ****************************************************************
void vtkSurface::ComputeSharpVertices(double treshold)
{

	// test first
	if (this->SharpVertices) return;

	double p1[3];
	double p2[3];
	double p3[3];
	double v1[3];
	double v2[3];
	double v3[3];
	double n;

	vtkIdType i;
	vtkIdType s1;
	vtkIdType s2;
	vtkIdType s3;
	vtkIdType f1;
	vtkIdType f2;
	vtkIdType j;
	bool ok;

	vtkPoints *points = this->GetPoints();

	vtkDoubleArray *FacesNormals = vtkDoubleArray::New();
	FacesNormals->SetNumberOfComponents(3);
	FacesNormals->SetNumberOfTuples(this->GetNumberOfCells());

	for (i=0;i<this->GetNumberOfCells();i++)
	{
		this->GetFaceVertices(i,s1,s2,s3);
		points->GetPoint(s1,p1);
		points->GetPoint(s2,p2);
		points->GetPoint(s3,p3);

		for (j=0;j<3;j++)
		{
			v1[j] = p2[j]-p1[j];
			v2[j] = p3[j]-p1[j];
		}

		v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
		v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
		v3[2] = v1[0]*v2[1] - v1[1]*v2[0];

		n = 1.0 / sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);

		v3[0] *= n;
		v3[1] *= n;
		v3[2] *= n;

		FacesNormals->SetTuple(i,v3);

	}

	this->SharpVertices = vtkIntArray::New();
	this->SharpVertices->SetNumberOfValues(this->GetNumberOfPoints());

	vtkIdList *vlist = vtkIdList::New();

	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		ok = 0;
		this->GetVertexNeighbourEdges(i,vlist);
		for (j=0;j<vlist->GetNumberOfIds();j++)
		{
			this->GetEdgeFaces(vlist->GetId(j),f1,f2);
			if (f2 >= 0)
			{
				FacesNormals->GetTuple(f1,v1);
				FacesNormals->GetTuple(f2,v2);

				v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
				v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
				v3[2] = v1[0]*v2[1] - v1[1]*v2[0];

				n = sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);


				if (n>treshold)
				{
					j  = 10000;
					ok = 1;
				}

			}
			else
			{
				j  = 10000;
				ok = 1;
			}
		}
		this->SharpVertices->SetValue(i,ok);
	}

	if (this->GetPointData()->GetScalars() != NULL) this -> GetPointData()->GetScalars()->Delete();
	this->GetPointData()->SetScalars(this->SharpVertices);

	FacesNormals->Delete();
	vlist->Delete();
}
// ****************************************************************

void vtkSurface::DeleteSharpVertices()
{
	if (this->SharpVertices)
	{
		this->SharpVertices->Delete();
		this->SharpVertices=0;
		this->GetPointData()->SetScalars(0);
	}
};


void vtkSurface::DeleteVerticesAreas()
{
	if (this->VerticesAreas)
	{
		this->VerticesAreas->Delete();
		this->VerticesAreas=0;
	}
};
void vtkSurface::DeleteEdgeLengths()
{
	if (this->EdgeLengths)
	{
		this->EdgeLengths->Delete();
		this->EdgeLengths=0;
	}
};
void vtkSurface::DeleteTrianglesAreas()
{
	if (this->TrianglesAreas)
	{
		this->TrianglesAreas->Delete();
		this->TrianglesAreas=0;
	}
};
void vtkSurface::DeleteTrianglesNormals()
{
	if (this->TrianglesNormals)
	{
		this->TrianglesNormals->Delete();
		this->TrianglesNormals=0;
	}
};


// **
// ** Le new de vtk
// **
// **
vtkSurface* vtkSurface::New()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkSurface");
	if(ret)
	{
		return (vtkSurface*)ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkSurface);
}

void vtkSurface::UnQuantizeCoordinates()
{
	vtkIdType i;
	double p[3];

	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->Points->GetPoint(i,p);

		p[0] = (p[0] / this->Factor) + this->Tx;
		p[1] = (p[1] / this->Factor) + this->Ty;
		p[2] = (p[2] / this->Factor) + this->Tz;

		this->Points->SetPoint(i,p);
	}
}

// ****************************************************************
// ****************************************************************
// fonction QuantizeCoordinates
// this function quantizes the vertices coordinates to integers coordinates
void vtkSurface::QuantizeCoordinates(int q)
{
	double xmin,ymin,zmin;
	double xmax,ymax,zmax;
	double p[3];
	double xdyn,ydyn,zdyn;
	double dyn;
	double factor;
	double Bounds[6];

	vtkIdType i;


	dyn = 0;

	this->GetPoints()->ComputeBounds();
	this->GetPoints()->GetBounds(Bounds);
	xmin=Bounds[0];
	xmax=Bounds[1];
	ymin=Bounds[2];
	ymax=Bounds[3];
	zmin=Bounds[4];
	zmax=Bounds[5];

	cout<<endl;

	xdyn = xmax-xmin;
	ydyn = ymax-ymin;
	zdyn = zmax-zmin;

	if (dyn<xdyn) dyn = xdyn;
	if (dyn<ydyn) dyn = ydyn;
	if (dyn<zdyn) dyn = zdyn;

	factor = pow(2.0,q) - 1;
	factor /= dyn;


	this->Tx=xmin;
	this->Ty=ymin;
	this->Tz=zmin;
	this->Factor=factor;
	/*
	this->QuantizeScaling(  this->Tx,this->Tx1,this->Tx2);
	this->UnQuantizeScaling(this->Tx,this->Tx1,this->Tx2);
	this->QuantizeScaling(  this->Ty,this->Ty1,this->Ty2);
	this->UnQuantizeScaling(this->Ty,this->Ty1,this->Ty2);

	this->QuantizeScaling(  this->Tz,this->Tz1,this->Tz2);
	this->UnQuantizeScaling(this->Tz,this->Tz1,this->Tz2);

	this->QuantizeScaling(  this->Factor,this->Factor1,this->Factor2);
	this->UnQuantizeScaling(this->Factor,this->Factor1,this->Factor2);

	factor = this->Factor;
	xmin   = this->Tx;
	ymin   = this->Ty;
	zmin   = this->Tz;
	*/

	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->Points->GetPoint(i,p);
		p[0] = floor( (p[0]-xmin) * factor + 0.5 );
		p[1] = floor( (p[1]-ymin) * factor + 0.5 );
		p[2] = floor( (p[2]-zmin) * factor + 0.5 );
		this->Points->SetPoint(i,p);
	}
	this->Modified();
}
// ****************************************************************
// ****************************************************************
// fonction QuantizeCoordinates
// this function quantizes the vertices coordinates to integers coordinates
void vtkSurface::QuantizeCoordinatesLike(vtkSurface *Mesh)
{
	int i;
	double p[3];
	Mesh->GetScalingFactors(this->Factor,this->Tx,this->Ty,this->Tz);
	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->Points->GetPoint(i,p);
		p[0] = floor( (p[0]-Tx) * Factor + 0.5 );
		p[1] = floor( (p[1]-Ty) * Factor + 0.5 );
		p[2] = floor( (p[2]-Tz) * Factor + 0.5 );
		this->Points->SetPoint(i,p);
	}
}

void vtkSurface::QuantizeCoordinates(double Factor, double Tx, double Ty, double Tz)
{
	this->Factor = Factor;
	this->Tx = Tx;
	this->Ty = Ty;
	this->Tz = Tz;

	int i;
	double p[3];
	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->Points->GetPoint(i,p);
		p[0] = floor( (p[0]-Tx) * Factor + 0.5 );
		p[1] = floor( (p[1]-Ty) * Factor + 0.5 );
		p[2] = floor( (p[2]-Tz) * Factor + 0.5 );
		this->Points->SetPoint(i,p);
	}
}

void vtkSurface::WriteToFile (const char *FileName)
{
	char filename[500];

	strcpy(filename,FileName);
	if (filename != NULL) {
		char *p;
		for (p = filename; *p; ++p)
			*p = tolower(*p);
	}

	if (strstr(filename,".vtk") != NULL) {
		vtkPolyDataWriter *Writer = vtkPolyDataWriter::New();
		Writer->SetInputData(this);
#if ( ( VTK_MAJOR_VERSION >= 9 ) && ( VTK_MINOR_VERSION >= 1 ) )
		Writer->SetFileVersion(42);
#endif
		Writer->SetFileName(FileName);
		Writer->Write();
	}

	if (strstr(filename,".ply") != NULL) {
		vtkPLYWriter *Writer = vtkPLYWriter::New();
		Writer->SetInputData(this);
		Writer->SetFileName(FileName);
		Writer->Write();
	}
}

// ****************************************************************
// ****************************************************************
void vtkSurface::CreateFromFile(const char *FileName)
{

	int ch;
	char *terminaison;
	char fin[180];
	char filename[180];

	strcpy(filename,FileName);

	if (filename != NULL)
	{
		char *p;

		for (p = filename; *p; ++p)
			*p = tolower(*p);
	}


	ch='.';
	strcpy (fin,".vtk");
	terminaison = strstr(filename,fin);
	if (terminaison!=NULL)
	{
		vtkPolyDataReader *Reader = vtkPolyDataReader::New();
		Reader->SetFileName(FileName);
		Reader->Update();
		this->CreateFromPolyData(Reader->GetOutput());
		Reader->Delete();
		return;
	}

	ch='.';
	strcpy (fin,".vtp");
	terminaison = strstr(filename,fin);
	if (terminaison!=NULL)
	{
		vtkPolyDataReader *Reader = vtkPolyDataReader::New();
		Reader->SetFileName(FileName);
		Reader->Update();
		this->CreateFromPolyData(Reader->GetOutput());
		Reader->Delete();
		return;
	}

	ch='.';
	strcpy (fin,".ply");
	terminaison = strstr(filename,fin);
	if (terminaison!=NULL)
	{
		vtkPLYReader *Reader = vtkPLYReader::New();
		Reader->SetFileName(FileName);
		Reader->Update();
		this->CreateFromPolyData(Reader->GetOutput());
		Reader->Delete();
		return;
	}

	ch='.';
	strcpy (fin,".wrl");
	terminaison=strstr(filename,fin);
	if (terminaison!=NULL)
	{
		cout<<"wrl"<<endl;

		vtkRenderer *vrmlrenderer1=vtkRenderer::New();
		vtkRenderWindow *vrmlrenWin=vtkRenderWindow::New();
		vrmlrenWin->AddRenderer(vrmlrenderer1);
		vtkVRMLImporter *importer=vtkVRMLImporter::New();
		importer->SetRenderWindow(vrmlrenWin);
		importer->SetFileName(FileName);
		importer->Read();
		vtkRendererCollection *renCollection=vrmlrenWin->GetRenderers();
		renCollection->InitTraversal();
		vtkRenderer *ren=renCollection->GetNextItem();
		vtkPolyData *temporary=(vtkPolyData*) ren->GetActors()->GetLastActor()->GetMapper()->GetInput();
		this->CreateFromPolyData(temporary);
		vrmlrenWin->Delete();
		vrmlrenderer1->Delete();
		importer->Delete();
	}
	ch='.';
	strcpy (fin,".off");
	terminaison = strstr(filename,fin);
	if (terminaison!=NULL)
	{
		vtkOFFReader *Reader = vtkOFFReader::New();
		Reader->SetFileName(FileName);
//		Reader->SetInfoOnCellsOn();
		Reader->Update();
		this->CreateFromPolyData(Reader->GetOutput());
		Reader->Delete();
		return;
	}
	ch='.';
	strcpy (fin,".obj");
	terminaison = strstr(filename,fin);
	if (terminaison!=NULL)
	{
		vtkOBJReader *Reader = vtkOBJReader::New();
		Reader->SetFileName(FileName);
		Reader->Update();
		this->CreateFromPolyData(Reader->GetOutput());
		Reader->Delete();
		return;
	}
	ch='.';
	strcpy (fin,".stl");
	terminaison = strstr(filename,fin);
	if (terminaison!=NULL)
	{
		vtkSTLReader *Reader = vtkSTLReader::New();
		Reader->SetFileName(FileName);
		Reader->Update();
		this->CreateFromPolyData(Reader->GetOutput());
		Reader->Delete();
		return;
	}

}

// ****************************************************************
// ****************************************************************
vtkSurface::vtkSurface()
:vtkSurfaceBase()
{
	this->SharpVertices=0;
	this->VerticesAreas=0;
	this->TrianglesAreas=0;
	this->TrianglesNormals=0;

	this->EdgeLengths=0;
	this->ConnectedComponents=0;
}

// ****************************************************************
// ****************************************************************
vtkSurface::~vtkSurface()
{
	this->DeleteSharpVertices();
	this->DeleteVerticesAreas();
	this->DeleteTrianglesAreas();
	this->DeleteTrianglesNormals();
	this->DeleteEdgeLengths();
	this->DeleteConnectedComponents();
}

// ****************************************************************
// ****************************************************************
