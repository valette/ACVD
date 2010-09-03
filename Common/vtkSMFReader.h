/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSMFReader.h,v $
  Language:  C++

  Copyright (c) 2003 Arnaud Gelas, Min-Su KIM 
  All rights reserved.

=========================================================================*/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Arnaud Gelas, Min-Su KIM
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

// .NAME vtkSMFReader - read Qsilm .smf files
// .SECTION Description
// vtkSMFReader is a source object that reads Qsilm .smf
// files. The output of this source object is polygonal data.
// .SECTION See Also
// vtkSMFImporter

#ifndef __vtkSMFReader_h
#define __vtkSMFReader_h

#include <vtkPolyDataSource.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>

class VTK_EXPORT vtkSMFReader : public vtkPolyDataSource 
{
public:
  static vtkSMFReader *New();
  
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of Qsilm .smf file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  vtkSMFReader();
  ~vtkSMFReader();
  
  void Execute();

  char *FileName;
private:
  vtkSMFReader(const vtkSMFReader&);  // Not implemented.
  void operator=(const vtkSMFReader&);  // Not implemented.

	int AddNormalVector(char *line, vtkFloatArray *NormalVector);
	int AddTextureCoordinate(char *line, vtkFloatArray *TextureTuple);
	int AddColorComponent(char *line,vtkFloatArray *colorTuple);
	int AddFace(char *line, vtkCellArray *polys);
	int AddVertex(char *line, vtkPoints *points);
};

#endif


