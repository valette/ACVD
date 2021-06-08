/*=========================================================================

  Program:   Mailleur 3D multi-r�solution (Creatis 2000 ~ nowadays)
  Module:    vtkOFFWriter.cxx
  Language:  C++
  Date:      2003/05
  Auteurs:   Alexandre Gouaillard
=========================================================================*/
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>

#include "vtkOFFWriter.h"
#include "vtkSurface.h"

//vtkCxxRevisionMacro(vtkOFFWriter, "$Revision: 1.5 $");
vtkStandardNewMacro(vtkOFFWriter);

//--------------------------------------------------------------------------
void vtkOFFWriter::WriteData()
{
	
	
	// Check input first
	if ( this->FileName == NULL)
    {
		vtkErrorMacro(<< "Please specify FileName to write");
		return;
    }

	FILE *fp;

    if ((fp = fopen(this->FileName, "w")) == NULL)
    {
		vtkErrorMacro(<< "Couldn't open file: " << this->FileName);
		return;
    }

	// if input is correct , proceed.
	int i;
	double P[3];

	vtkPolyData *polydata = this->GetInput();
	vtkSurface *surfaceTemp = vtkSurface::New();
	surfaceTemp->CreateFromPolyData(polydata);

	int nbpts = polydata->GetNumberOfPoints();
	int nbcls = polydata->GetNumberOfCells();
	int nbEdg = surfaceTemp->GetNumberOfEdges();

	fprintf (fp,"OFF\n");
	fprintf (fp,"%d %d %d\n", nbpts, nbcls,nbEdg);

	// * write out the coordinates *
	for (i = 0; i < nbpts; i++)
	{
		polydata->GetPoint(i,P);
		fprintf (fp,"%f %f %f\n", P[0],P[1],P[2]);
	}

	vtkIdType nbPtsCell;
#if ( (VTK_MAJOR_VERSION < 9))
		vtkIdType *ptIdList;
#else
		const vtkIdType *ptIdList;
# endif
	vtkIdType j;

	// * write the faces *
	for (i = 0; i <nbcls;i++)
	{
		polydata->GetCellPoints(i,nbPtsCell,ptIdList);
		int temp = nbPtsCell;
		fprintf (fp,"%d ", temp);
		
		for (j=0; j<nbPtsCell ;j++)
		{
			temp = ptIdList[j];
			fprintf (fp,"%d ", temp);
		}
		
		fprintf (fp,"\n");
	
	}

	fclose(fp);

	surfaceTemp->Delete();
}
//--------------------------------------------------------------------------


//----------------------------------------------------------------------------
/*void vtkOFFWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}*/
