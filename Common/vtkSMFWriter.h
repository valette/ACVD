/*=========================================================================

  Program:   SMF Writer
  Module:    vtkSMFWriter.cxx
  Language:  C++
  Date:      2003/05
  Auteurs:   Sebastien Valette
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

// .NAME vtkSMFWriter - write SMF files
// .SECTION Description
// vtkSMFWriter writes SMF files (.smf) files.

#ifndef __vtkSMFWriter_h
#define __vtkSMFWriter_h

#include <vtkPolyDataWriter.h>

class VTK_EXPORT vtkSMFWriter : public vtkPolyDataWriter
{
public:
	static vtkSMFWriter *New();
	vtkTypeMacro(vtkSMFWriter,vtkPolyDataWriter);
	//virtual void PrintSelf(ostream& os, vtkIndent indent);

protected:
	vtkSMFWriter() {};
	~vtkSMFWriter() {};

	void WriteData()
	{
		std::ofstream File;
		File.open (this->FileName, std::ofstream::out | std::ofstream::trunc);

		vtkIdType i;
		vtkIdType nverts=this->GetInput()->GetNumberOfPoints();
		vtkIdType nfaces=this->GetInput()->GetNumberOfCells();

		File << "# Generated from vtkSMFWriter" << endl;
		File<< "# " <<  nverts<< " vertices" << endl;
		File << "# " << nfaces << " faces" << endl;

		double P[3];
		for(i=0; i<nverts; i++)
		{
			this->GetInput()->GetPoint(i,P);
			File << "v "<< P[0] << " " << P[1] << " " << P[2] << endl;
		}

		const vtkIdType *Vertices;
		vtkIdType j,NumberOfVertices;
		for(i=0; i<nfaces; i++)
		{
			this->GetInput()->GetCellPoints(i,NumberOfVertices,Vertices);
			File << "f ";
			for (j=0;j<NumberOfVertices;j++)
				File<< Vertices[j]+1 << " ";
			File<< endl;
		}
		File.close();
	}

	private:
	vtkSMFWriter(const vtkSMFWriter&);  // Not implemented.
	void operator=(const vtkSMFWriter&);  // Not implemented.
};
vtkStandardNewMacro(vtkSMFWriter);
#endif

