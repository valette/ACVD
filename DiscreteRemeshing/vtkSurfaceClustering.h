/***************************************************************************
                          vtkSurfaceClustering.h  -  description
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
#ifndef _VTKSURFACECLUSTERING_H_
#define _VTKSURFACECLUSTERING_H_

#include <vtkTriangleFilter.h>

#include "vtkUniformClustering.h"
#include "vtkThreadedClustering.h"
#include "vtkLloydClustering.h"

/**
 * This class provides an abstract layer to interface vtkUniformClustering with vtkSurface objects
 * Two classes derive from vtkSurfaceClustering : vtkVerticesProcessing and vtkTrianglesProcessing
 */

#ifdef DOmultithread
template < class Metric > class vtkSurfaceClustering:public vtkThreadedClustering <Metric>
#else
#ifdef DOLloydClustering
template < class Metric > class vtkSurfaceClustering:public vtkLloydClustering <Metric>
#else
template < class Metric > class vtkSurfaceClustering:public vtkUniformClustering <Metric>
#endif
#endif
{

public:

	/// Sets the Input mesh
	void SetInput(vtkSurface *Input);
	
	/// Returns the Input mesh
	vtkGetMacro( Input, vtkSurface* );

	/// Returns the Rendering window used for display
	RenderWindow* GetDisplayWindow(){ return this->Window; };

	vtkIdType GetNumberOfEdges() { return this->Input->GetNumberOfEdges(); }

	void GetEdgeItemsSure (vtkIdType Item, vtkIdList * VList) {

		if ( this->ClusteringType == 0 )
			this->GetInput()->GetEdgeFaces( Item, VList );
		else {

			vtkIdType v1, v2;
			VList->Reset();
			this->Input->GetEdgeVertices ( Item, v1, v2 );
			VList->InsertNextId( v1 );
			VList->InsertNextId( v2 );

		}

	};

	vtkGetMacro( ClusteringType, int );

protected:

	/// Parameter indicating the type of Clustering (0: Faces ;1:Vertices)
	int ClusteringType;
	virtual int GetNumberOfDualItems() = 0;
	void CreateWindows();
	void InteractWithClusteringWindow();
	void Snapshot();
	void BuildMetric();

	vtkSurfaceClustering ();
	~vtkSurfaceClustering ();
	vtkSurface *Input;	// input vtkSurface
	RenderWindow *Window; // The window where the clustering is displayed

};



template <class Metric>
void vtkSurfaceClustering<Metric>::BuildMetric()
{

	this->MetricContext.BuildMetric( this->Input, this->ClusteringType );

}

template <class Metric>
void vtkSurfaceClustering<Metric>::CreateWindows()
{
	if (this->Display)
	{
		Window=RenderWindow::New();
		Window->SetInputData(this->Input);
		if (this->AnchorRenderWindow)
			Window->AttachToRenderWindow(this->AnchorRenderWindow);
		else
			this->AnchorRenderWindow=Window;

		if (this->ClusteringType==1)
		{
			this->Input->GetCellData()->SetScalars(0);
			this->Input->GetPointData()->SetScalars(this->Clustering);
		}
		else
		{
			this->Input->GetCellData()->SetScalars(this->Clustering);
			this->Input->GetPointData()->SetScalars(0);
		}

		Window->DisplayRandomColors(this->NumberOfClusters+1);
		Window->Render();
		Window->SetWindowName("Clustering");

		Window->Interact();


		if (this->Display>2)
		{
			std::stringstream strfile;
			strfile<<"movie/mesh"<<(this->FrameNumber++)+1000<<".png";
			Window->Capture(strfile.str().c_str());
		}
	}
}

template <class Metric>
void vtkSurfaceClustering<Metric>::Snapshot()
{
	if (this->Display)
	{
		this->Clustering->Modified();
		Window->Render();
		if (this->Display>2)
		{
			std::stringstream strfile;
			strfile<<"movie/mesh"<<(this->FrameNumber++)+1000<<".png";
			Window->Capture(strfile.str().c_str());
		}		
	}
}
template <class Metric>
void vtkSurfaceClustering<Metric>::InteractWithClusteringWindow()
{
	if (this->Display)
	{
		this->Window->Render();
		this->Window->Interact();
	}
}


template <class Metric>
void vtkSurfaceClustering<Metric>::SetInput(vtkSurface *Input)
{
	if (this->Input)
		this->Input->UnRegister(this);

	vtkIdType NV;
	vtkIdType *Vertices;
	bool trianglesOnly = true;
	for (vtkIdType i = 0; i < Input->GetNumberOfCells(); i++) {
		if (Input->IsFaceActive(i)) {
			Input->GetFaceVertices(i, NV, Vertices);
			if (NV != 3) {
				trianglesOnly = false;
				break;
			}
		}
	}

	if (!trianglesOnly) {
		// triangulate the mesh if it contains polygons
		vtkTriangleFilter *triangulate = vtkTriangleFilter::New();
		triangulate->SetInputData(Input);
		triangulate->PassVertsOff ();
		triangulate->Update();
		Input = vtkSurface::New();
		Input->CreateFromPolyData(triangulate->GetOutput());
		triangulate->Delete();
	}

	this->Input=Input;
	if (Input)
		this->Input->Register(this);
};


template <class Metric> vtkSurfaceClustering<Metric>::vtkSurfaceClustering ()
{
	this->AnchorRenderWindow=0;
	this->Input=0;
	this->Window=0;
}

template <class Metric> vtkSurfaceClustering<Metric>::~vtkSurfaceClustering ()
{
	if (this->Input)
		this->Input->Delete();
		
	if (this->Window)
		this->Window->Delete();
}
#endif
