/*=========================================================================

  Program:   Visualization Toolkit
  Language:  C++
  Thanks:    Thanks to Abdalmajeid M. Alyassin who developed this class.


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

#include <vtkObjectFactory.h>
#include <vtkCommand.h>
#include <vtkMath.h>
#include <vtkTriangle.h>
#include <vtkIdList.h>
#include <vtkCell.h>


#include "vtkVolumeProperties.h"

//------------------------------------------------------------------------------
vtkVolumeProperties* vtkVolumeProperties::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkVolumeProperties");
  if(ret)
    {
    return (vtkVolumeProperties*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkVolumeProperties;
}




#define  VTK_CUBE_ROOT(x) ((x<0.0)?(-pow((-x),0.333333333333333)):(pow((x),0.333333333333333)))

// Constructs with initial 0 values.
vtkVolumeProperties::vtkVolumeProperties()
{
  this->SurfaceArea = 0.0;
  this->Volume  = 0.0;
  this->SignedVolume  = 0.0;
  this->XG = 0.0;
  this->YG = 0.0;
  this->ZG = 0.0;
  this->NormalizedShapeIndex = 0.0;
  this->InputData = 0;
}

// Constructs with initial 0 values.
vtkVolumeProperties::~vtkVolumeProperties()
{
}


// Description:
// Make sure input is available then call up execute method...
void vtkVolumeProperties::Update()
{
  vtkPolyData *input = this->GetInputData();
  
  // make sure input is available
  if ( ! input )
    {
    vtkErrorMacro(<< "No input...can't execute!");
    return;
    }

  if (input->GetMTime() > this->ExecuteTime || 
      this->GetMTime() > this->ExecuteTime ) {
    this->InvokeEvent(vtkCommand::StartEvent,NULL);

    this->Execute();
    this->ExecuteTime.Modified();
    this->InvokeEvent(vtkCommand::EndEvent,NULL);
    }
}

// Description:
// This method measures volume, surface area, and normalized shape index.
// Currently, the input is a ploydata which consists of triangles.
void vtkVolumeProperties::Execute()
{
	vtkIdList *ptIds;
	vtkPolyData *input = this->GetInputData();
	int cellId, numCells, numPts, numIds;
	double *p;

	numCells=input->GetNumberOfCells();
	numPts = input->GetNumberOfPoints();
	if (numCells < 1 || numPts < 1)
	{
		vtkErrorMacro(<<"No data to measure...!");
		return;
	}

	ptIds = vtkIdList::New();
	ptIds->Allocate(VTK_CELL_SIZE);

	//
	// Traverse all cells, obtaining node coordinates.
	//
	double    u[3],v[3],w[3],bary[3],m[3];
	double    P1[3],P2[3],P3[3];
	int      idx;

	for ( idx =0; idx < 3 ; idx++ ) 
		bary[idx]  = 0.0;

	//coordonnées du barycentre du maillage
	for ( idx=0; idx < numPts ; idx++)
	{
		p = input->GetPoint(idx);
		bary[0] += p[0];
		bary[1] += p[1];
		bary[2] += p[2];
	}

	bary[0] /= (double) numPts ;
	bary[1] /= (double) numPts ;
	bary[2] /= (double) numPts ;

	this->XG = bary[0] ;
	this->YG = bary[1] ;
	this->ZG = bary[2] ;
	
	this->SurfaceArea=0;
	this->SignedVolume=0;

	//parcours de toutes les cellules
	for (cellId=0; cellId < numCells; cellId++)
	{
		if ( input->GetCellType(cellId) != VTK_TRIANGLE)
		{
			vtkErrorMacro(<<"Sorry Data type has to be VTK_TRIANGLE only " << input->GetCellType(cellId));
			return;
		}
		input->GetCellPoints(cellId,ptIds);
		numIds = ptIds->GetNumberOfIds();
		input->GetPoint(ptIds->GetId(0),P1);
		input->GetPoint(ptIds->GetId(1),P2);
		input->GetPoint(ptIds->GetId(2),P3);
		
		for (int i=0;i<3;i++)
		{
			u[i]=P1[i]-bary[i];
			v[i]=P2[i]-bary[i];
			w[i]=P3[i]-bary[i];
		}
		
		vtkMath::Cross(v,w,m);
		this->SignedVolume+=vtkMath::Dot(m,u)/6.0;
		this->SurfaceArea+= vtkTriangle::TriangleArea(P1,P2,P3);
		
	}

	this->Volume = fabs(this->SignedVolume) ;
	this->NormalizedShapeIndex = (sqrt(this->SurfaceArea)/VTK_CUBE_ROOT(this->Volume))/2.199085233;
	ptIds->Delete();
}

/*void vtkVolumeProperties::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkProcessObject::PrintSelf(os,indent);

  if (!this->GetInput()) 
    {
    return;
    }
  os << indent << "VolumeC: " << this->GetVolumeC () << "\n";
  os << indent << "VolumeB: " << this->GetVolumeB () << "\n";
  os << indent << "XG: " << this->GetXG () << "\n";
  os << indent << "YG: " << this->GetYG () << "\n";
  os << indent << "ZG: " << this->GetZG () << "\n";
  os << indent << "Surface Area: " << this->GetSurfaceArea () << "\n";
  os << indent << "Normalized Shape Index C: " << this->GetNormalizedShapeIndexC () << "\n";
  os << indent << "Normalized Shape Index B: " << this->GetNormalizedShapeIndexB () << "\n";
}*/
