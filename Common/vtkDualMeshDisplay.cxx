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

#include <vtkCellData.h>
#include <vtkObjectFactory.h>
#include "vtkDualMeshDisplay.h"

RenderWindow* vtkDualMeshDisplay::GetRenderWindow()
{
	return (this->Window);
}

void vtkDualMeshDisplay::SetInputData(vtkSurface *Input)
{
	if (this->Input)
	{
		this->Input->UnRegister(this);
		if (this->Dual)
		{
			Dual->Delete();
			Dual=0;
		}
	}
	
	this->Input=Input;
	this->Input->Register(this);
}

void vtkDualMeshDisplay::SetColors(vtkIntArray *Colors)
{
	if (this->Colors)
		this->Colors->UnRegister(this);
	
	this->Colors=Colors;
	this->Colors->Register(this);
}

vtkSurface *vtkDualMeshDisplay::GetOutput()
{
	this->Update();
	return (this->Dual);
}

void vtkDualMeshDisplay::Update()
{
	if (this->Dual==0)
		this->CreateDual();

	if (this->Colors==0)
		return;

	int NumFaces=this->Input->GetNumberOfCells();
	vtkIntArray *Colors=(vtkIntArray*) this->Dual->GetCellData()->GetScalars();
	for (int i=0;i<NumFaces;i++)
	{
		vtkIdType v1,v2,v3;
		this->Input->GetFaceVertices(i,v1,v2,v3);
		Colors->SetValue(i*6,this->Colors->GetValue(v1));
		Colors->SetValue(i*6+1,this->Colors->GetValue(v1));
		Colors->SetValue(i*6+2,this->Colors->GetValue(v2));
		Colors->SetValue(i*6+3,this->Colors->GetValue(v2));
		Colors->SetValue(i*6+4,this->Colors->GetValue(v3));
		Colors->SetValue(i*6+5,this->Colors->GetValue(v3));
	}
	Colors->Modified();
}

void vtkDualMeshDisplay::CreateDual()
{
	int NumFaces=this->Input->GetNumberOfCells();
	int NumPoints=this->Input->GetNumberOfPoints();
	int NumEdges=this->Input->GetNumberOfEdges();
	vtkIdType i;
	vtkIdType v1,v2,v3;

	double Point1[3];
	double Point2[3];
	double Point3[3];

	this->Dual=vtkSurface::New();
	this->Window->SetInputData(this->Dual);

	for( i = 0; i < NumPoints; i++ )
	{
		this->Input->GetPoint( i, Point1 );
		Dual->AddVertex (Point1 );
	}

	for( i = 0; i < NumEdges; i++ )
	{
		if( this->Input->IsEdgeActive( i ) )
		{
			this->Input->GetEdgeVertices( i, v1, v2 );
			this->Input->GetPoint( v1, Point1 );
			this->Input->GetPoint( v2, Point2 );
			for( int j = 0; j < 3; j++ )
				Point1[j] = 0.5 * ( Point1[j] + Point2[j] );

			this->Dual->AddVertex(Point1);
		}
	}

	for( i = 0; i < NumFaces; i++ )
	{
		if( this->Input->IsFaceActive( i ) )
		{
			this->Input->GetFaceVertices( i, v1, v2, v3 );
			this->Input->GetPoint(v1,Point1);
			this->Input->GetPoint(v2,Point2);
			this->Input->GetPoint(v3,Point3);
			for (int j=0;j<3;j++)
				Point1[j]=(Point1[j]+Point2[j]+Point3[j])/3.0;
			this->Dual->AddVertex(Point1);
		}
	}

	for (i=0;i<NumFaces;i++)
	{
		this->Input->GetFaceVertices(i,v1,v2,v3);
		vtkIdType v4,v5,v6,v7;
		v4=this->Input->IsEdge(v1,v2)+NumPoints;
		v5=this->Input->IsEdge(v2,v3)+NumPoints;
		v6=this->Input->IsEdge(v3,v1)+NumPoints;
		v7=i+NumPoints+NumEdges;
		this->Dual->AddFace(v1,v4,v7);
		this->Dual->AddFace(v1,v7,v6);
		this->Dual->AddFace(v2,v5,v7);
		this->Dual->AddFace(v2,v7,v4);
		this->Dual->AddFace(v3,v7,v5);
		this->Dual->AddFace(v3,v6,v7);
	}

	this->Dual->Modified();
	if (this->Colors)
	{
		vtkIntArray *Array=vtkIntArray::New();
		Array->SetNumberOfValues(NumFaces*6);
		this->Dual->GetCellData()->SetScalars(Array);
		Array->Delete();
	}
}

vtkDualMeshDisplay *
vtkDualMeshDisplay::New ()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject *ret = vtkObjectFactory::CreateInstance ("vtkDualMeshDisplay");
	if (ret)
	{
		return (vtkDualMeshDisplay *) ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkDualMeshDisplay);
}

vtkDualMeshDisplay::~vtkDualMeshDisplay()
{
	if (this->Input)
		this->Input->UnRegister(this);

	Window->Delete();
}

vtkDualMeshDisplay::vtkDualMeshDisplay()
{
	this->Input=0;
	this->Colors=0;
	this->Dual=0;
	
	this->Window=RenderWindow::New();
}
