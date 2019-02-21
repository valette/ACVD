/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSMFReader.cxx,v $
  Language:  C++

  Copyright (c) 2003 Arnaud Gelas, Min-Su KIM 
  All rights reserved.

=========================================================================*/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Authors : Arnaud Gelas, Min-Si KIM
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

#include "vtkSMFReader.h"

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include "vtkCellData.h"


// Description:
// Instantiate object with NULL filename.
vtkSMFReader* vtkSMFReader::New()
{
	  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkSMFReader");
  if(ret) return (vtkSMFReader*)ret;

  // If the factory was unable to create the object, then create it here.
  return new vtkSMFReader;
}

vtkSMFReader::vtkSMFReader()
{
  this->FileName = NULL;
}

vtkSMFReader::~vtkSMFReader()
{
  if (this->FileName)
    {
    delete [] this->FileName;
    this->FileName = NULL;
    }
}

int vtkSMFReader::AddVertex(char *line, vtkPoints *points)
{
	// this is a vertex definition, expect three floats, separated by whitespace:
	float xyz[3];

	int everything_ok = 1;

	if (sscanf(line, "v %f %f %f", xyz, xyz + 1, xyz + 2)==3) 
	{
		points->InsertNextPoint(xyz);
	}
	else 
	{
		vtkErrorMacro(<<"Error in reading file");
		everything_ok=0; // (false)
	}
	return everything_ok;
}

int vtkSMFReader::AddFace(char *line, vtkCellArray *polys)
{
	// this is a face definition, expect three integer, separated by whitespace:
	int temp[3];
	vtkIdList *iFace = vtkIdList::New();
	iFace->SetNumberOfIds(3);
	
	int everything_ok = 1;

	if (sscanf(line, "f %d %d %d", temp, temp+1, temp+2)==3) 
	{
		for(int i=0; i<3; i++)
			iFace->SetId(i,temp[i]-1);
		polys->InsertNextCell(iFace);
	}
	else 
	{
		vtkErrorMacro(<<"Error in reading file");
		everything_ok=0; // (false)
	}
	return everything_ok;
}

int vtkSMFReader::AddColorComponent(char *line,vtkFloatArray *colorTuple)
{
	// this is a face definition, expect three integer, separated by whitespace:
	float Color[3];

	int everything_ok = 1;

	if (sscanf(line, "c %f %f %f", Color, Color + 1, Color + 2)==3) 
	{
		for (int i=0; i<3 ; i++)
		{
			if((Color[i] > 1)||(Color[i] < 0))
			{
				vtkErrorMacro(<<"Error in reading file");
				everything_ok=0; // (false)
				break;
			}
			else
			{
				colorTuple->InsertNextTuple(Color);
			}
		}
	}
	else 
	{
		vtkErrorMacro(<<"Error in reading file");
		everything_ok=0; // (false)
	}
	
	return everything_ok;
}

int vtkSMFReader::AddTextureCoordinate(char *line, vtkFloatArray *TextureTuple)
{
	// this is a face definition, expect three integer, separated by whitespace:
	float Texture[3];
	
	int everything_ok = 1;

	if (sscanf(line, "c %f %f", Texture, Texture + 1)==2) 
	{
		for (int i=0; i<2 ; i++)
		{
			if((Texture[i] > 1)||(Texture[i] < 0))
			{
				vtkErrorMacro(<<"Error in reading file");
				everything_ok=0; // (false)
				break;
			}
			else
			{
				TextureTuple->InsertNextTuple(Texture);
			}
		}
	}
	else 
	{
		vtkErrorMacro(<<"Error in reading file");
		everything_ok=0; // (false)
	}
	return everything_ok;
}

int vtkSMFReader::AddNormalVector(char *line, vtkFloatArray *NormalVector)
{
	// this is a normal, expect three floats, separated by whitespace:
    float Normal[3];
	
	int everything_ok = 1;

	if (sscanf(line, "n %f %f %f", Normal, Normal + 1, Normal + 2)==3) 
    {
		float sum = 0.0;
		for (int i=0; i<3 ; i++)
		{
			sum += Normal[i] * Normal[i];
		}

		if(sum != 1)
		{
			vtkErrorMacro(<<"Error in reading file");
			everything_ok=0; // (false)
		}
		else
		{
			NormalVector->InsertNextTuple(Normal);
		}
    }
	else 
	{
		vtkErrorMacro(<<"Error in reading file");
		everything_ok=0; // (false)
	}

	return everything_ok;
}


void vtkSMFReader::Execute()
{
  if (!this->FileName) 
    {
    vtkErrorMacro(<< "A FileName must be specified.");
    return;
    }
    
  FILE *in = fopen(this->FileName,"r");
    
  if (in == NULL) 
    {
    vtkErrorMacro(<< "File " << this->FileName << " not found");
    return;
    }
    
  vtkDebugMacro(<<"Reading file");
    
  // intialise some structures to store the file contents in
  vtkPoints *points = vtkPoints::New(); 
  
  vtkFloatArray *tcoords = vtkFloatArray::New();
  tcoords->SetNumberOfComponents(2);
  
  vtkFloatArray *normals = vtkFloatArray::New();
  normals->SetNumberOfComponents(3);

  vtkCellArray *polys = vtkCellArray::New();

  vtkFloatArray *colorTuple = vtkFloatArray::New();
  colorTuple->SetNumberOfComponents(3);

  
  int everything_ok = 1; // (true)   (use of this flag avoids early return and associated memory leak)

  // -- work through the file line by line, assigning into the above six structures as appropriate --
  int FlagForVertex = 0;


	const int MAX_LINE=1024;
	char line[MAX_LINE];
	

	while (everything_ok && fgets(line,MAX_LINE,in)!=NULL) 
    {
		if ((strncmp(line,"bind c vertex ",14)==0) || (strncmp(line,"bind r vertex ",14)) || (strncmp(line,"bind n vertex ",14)) )
		{
			FlagForVertex = 1;
		}
	  
		if ((strncmp(line,"bind c face ",12)==0) || (strncmp(line,"bind r face ",12)) || (strncmp(line,"bind n face ",12)) )
		{
			FlagForVertex = 0;
		}

		// in the SMF format the first characters determine how to interpret the line:
		if (strncmp(line,"v",1)==0) 
		{
			everything_ok = this->AddVertex(line, points);
		}

		else if (strncmp(line,"f",1)==0)
		{
			everything_ok = this->AddFace(line, polys);
		}
		else if (strncmp(line,"r",1)==0)
		{
			everything_ok = this->AddTextureCoordinate(line, tcoords);
		}
		else if (strncmp(line,"n",1)==0)
		{
			everything_ok = this->AddNormalVector(line, normals);
		}
		else if (strncmp(line,"c",1)==0)
		{
			everything_ok = this->AddColorComponent(line, colorTuple);
		}

	}


    
  // we have finished with the file
  fclose(in); 

  if (everything_ok)   // (otherwise just release allocated memory and return)
    {
    // -- now turn this lot into a useable vtkPolyData --

    // if there are no tcoords or normals or they match exactly 
    // then we can just copy the data into the output (easy!)
    vtkDebugMacro(<<"Copying file data into the output directly");

     this->GetOutput()->SetPoints(points);
      this->GetOutput()->SetPolys(polys);

      // if there is an exact correspondence between tcoords and vertices then can simply
      // assign the tcoords points as point data

	  if(FlagForVertex == 1)
	  {
		  if (tcoords->GetNumberOfTuples() == points->GetNumberOfPoints())
		  {
			  this->GetOutput()->GetPointData()->SetTCoords(tcoords);
		  }
		  if (normals->GetNumberOfTuples() == points->GetNumberOfPoints())
		  {
			  this->GetOutput()->GetPointData()->SetNormals(normals);
		  }

	  }
	  else
	  {
		if (tcoords->GetNumberOfTuples() == polys->GetNumberOfCells())
		  {
			  this->GetOutput()->GetCellData()->SetTCoords(tcoords);
		  }
		  if (normals->GetNumberOfTuples() == polys->GetNumberOfCells())
		  {
			  this->GetOutput()->GetCellData()->SetNormals(normals);
		  }
	  }
      
      this->GetOutput()->Squeeze();
      
  }
    // otherwise we can duplicate the vertices as necessary (a bit slower)
   /* else 
      {
      vtkDebugMacro(<<"Duplicating vertices so that tcoords and normals are correct");

      vtkPoints *new_points = vtkPoints::New();
      vtkFloatArray *new_tcoords = vtkFloatArray::New();
      new_tcoords->SetNumberOfComponents(2);
      vtkFloatArray *new_normals = vtkFloatArray::New();
      new_normals->SetNumberOfComponents(3);
      vtkCellArray *new_polys = vtkCellArray::New();

      // for each poly, copy its vertices into new_points (and point at them)
      // also copy its tcoords into new_tcoords
      // also copy its normals into new_normals
      polys->InitTraversal();
      tcoord_polys->InitTraversal();
      normal_polys->InitTraversal();
      int i,j;
      vtkIdType dummy_warning_prevention_mechanism[1];
      vtkIdType n_pts=-1,*pts=dummy_warning_prevention_mechanism;
      vtkIdType n_tcoord_pts=-1,*tcoord_pts=dummy_warning_prevention_mechanism;
      vtkIdType n_normal_pts=-1,*normal_pts=dummy_warning_prevention_mechanism;
      for (i=0;i<polys->GetNumberOfCells();i++) 
        {
        polys->GetNextCell(n_pts,pts); 
        tcoord_polys->GetNextCell(n_tcoord_pts,tcoord_pts);
        normal_polys->GetNextCell(n_normal_pts,normal_pts);

        // If some vertices have tcoords and not others (likewise normals)
        // then we must do something else VTK will complain. (crash on render attempt)
        // Easiest solution is to delete polys that don't have complete tcoords (if there 
        // are any tcoords in the dataset) or normals (if there are any normals in the dataset).

        if ( (n_pts!=n_tcoord_pts && hasTCoords) || (n_pts!=n_normal_pts && hasNormals) ) 
          {
          // skip this poly
          vtkDebugMacro(<<"Skipping poly "<<i+1<<" (1-based index)");
          }
        else 
          {
          // copy the corresponding points, tcoords and normals across
          for (j=0;j<n_pts;j++) 
            {
            // copy the tcoord for this point across (if there is one)
            if (n_tcoord_pts>0)
              new_tcoords->InsertNextTuple(tcoords->GetTuple(tcoord_pts[j]));
            // copy the normal for this point across (if there is one)
            if (n_normal_pts>0)
              new_normals->InsertNextTuple(normals->GetTuple(normal_pts[j]));
            // copy the vertex into the new structure and update
            // the vertex index in the polys structure (pts is a pointer into it)
            pts[j] = new_points->InsertNextPoint(points->GetPoint(pts[j]));
            }
          // copy this poly (pointing at the new points) into the new polys list 
          new_polys->InsertNextCell(n_pts,pts);
          }
        }

      // use the new structures for the output
      this->GetOutput()->SetPoints(new_points);
      this->GetOutput()->SetPolys(new_polys);
      if (hasTCoords)
        this->GetOutput()->GetPointData()->SetTCoords(new_tcoords);
      if (hasNormals)
        this->GetOutput()->GetPointData()->SetNormals(new_normals);
      this->GetOutput()->Squeeze();

      new_points->Delete();
      new_polys->Delete();
      new_tcoords->Delete();
      new_normals->Delete();
      }
    }*/

  points->Delete();
  tcoords->Delete();
  normals->Delete();
  polys->Delete();
  //tcoord_polys->Delete();
  //normal_polys->Delete();
  colorTuple->Delete();
  


}

void vtkSMFReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
     << (this->FileName ? this->FileName : "(none)") << "\n";

}

