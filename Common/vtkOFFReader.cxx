/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOFFReader.cxx,v $
  Language:  C++
  Authors:   S. valette       -> basic code , for triangulations only
			 A. gouaillard	  -> debug mode, general code for any polygonal mesh

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

#include <vtkByteSwap.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkObjectFactory.h>

#include "vtkOFFReader.h"

#ifdef read
#undef read
#endif

#ifdef close
#undef close
#endif

#define ALEX_DEBUG 0

//----------------------------------------------------------------------------
vtkOFFReader* vtkOFFReader::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkOFFReader");
  if(ret) return (vtkOFFReader*)ret;

  // If the factory was unable to create the object, then create it here.
  return new vtkOFFReader;
}


//----------------------------------------------------------------------------
vtkOFFReader::vtkOFFReader()
{
  this->FileName = NULL;

  this->SetInfoOnCellsOff();
  this->NumberOfPoints = 0;
  this->NumberOfCells  = 0;
}

//----------------------------------------------------------------------------
vtkOFFReader::~vtkOFFReader()
{   
  if (this->FileName)
    {
    delete [] this->FileName;
    this->FileName = NULL;
    }
}

//----------------------------------------------------------------------------
// 
void vtkOFFReader::RequestInformation()
{
  vtkPolyData *output = this->GetOutput();

 // output->SetMaximumNumberOfPieces(-1);
}


//----------------------------------------------------------------------------
// 
void vtkOFFReader::RequestData()
{
	vtkPoints *temp_points=vtkPoints::New();
	vtkPolyData *output_mesh=this->GetOutput();
	vtkCellArray *temp_cells=vtkCellArray::New();

	FILE *stream;
	float fbuffer;
	float temp_point[3];
	vtkIdType i;
	vtkIdType j;
	int buffer;
	vtkIdType Temp_Cell[200];
	vtkIdType CellType;
	char format[81];

	// this variable is only used to avoid compilation warnings...
	int UnusedResult;
	
	float tampon; 

	// ouverture pas propre
	stream = fopen(this->FileName, "r" );
	if( stream == NULL ) {cerr << "The file " << this->FileName << "was not opened" << endl;}
	else
	 {
		// recuperation du type du fichier, et on jete :)
		UnusedResult=fscanf( stream, "%s", format );

		// recuperation de la ligne d info : nb de points
		UnusedResult=fscanf( stream, "%d", &buffer);
		this->NumberOfPoints = buffer;

#if ALEX_DEBUG 
		cout << "nb de points " << buffer << endl;
#endif

		// allocation memoire pour les points aue nous allons lire
		temp_points->Allocate(this->NumberOfPoints);

		// recuperation de la ligne d info : nb de cellules
		UnusedResult=fscanf( stream, "%d", &buffer);
		this->NumberOfCells=buffer;
#if ALEX_DEBUG
		cout << "nb of cellules " << buffer << endl;
#endif

		// allocation memoire pour les cellules que nous allons lire
		temp_cells->Allocate(this->NumberOfCells);

		// lecture de la ligne d info : nombre d'arretes
		UnusedResult=fscanf( stream, "%d", &buffer );

		// boucle sur les points : lecture des coordonnees
		for (i=0;i<this->NumberOfPoints;i++)
		{

#if ALEX_DEBUG			
			cout << "point en cours de chargement :" << i << endl;
#endif
			
			// coordonees en X
			UnusedResult=fscanf( stream, "%f", &fbuffer);
			temp_point[0]= fbuffer;
			
#if ALEX_DEBUG
			cout << "X = " << fbuffer;
#endif

			// coordonnees en Y
			UnusedResult=fscanf( stream, "%f", &fbuffer);
			temp_point[1]= fbuffer;

#if ALEX_DEBUG
			cout << "Y = " << fbuffer;
#endif

			// coordonnees en Z
			UnusedResult=fscanf( stream, "%f", &fbuffer);
			temp_point[2]= fbuffer;

#if ALEX_DEBUG
			cout << "Z = " << fbuffer;
#endif

			// creation du point dans notre liste
			temp_points->InsertPoint(i,temp_point);

#if ALEX_DEBUG
			cout << endl;
#endif
		}

		// boucle sur les cellules : lecture des coordonnees 
		for (i=0;i<this->NumberOfCells;i++)
		{
			int intValue;
			// lecture du nombre de points de la cellule
			UnusedResult=fscanf( stream, "%d", &intValue);
			CellType = intValue;

#if ALEX_DEBUG
			cout << "cellule a " << CellType << " points" << endl;
#endif
			
			// lecture des IDs des points de la cellule
			for (j=0;j<CellType;j++) 
			{
				UnusedResult=fscanf( stream, "%d", &intValue);
				Temp_Cell[j] = intValue;

#if ALEX_DEBUG
				cout << "point # " << j << " Id " << Temp_Cell[j] << "; "; 
#endif
			
			}
#if ALEX_DEBUG
			cout << endl;
#endif			

			// creation de la cellule dans notre liste
			temp_cells->InsertNextCell(CellType,Temp_Cell);
			
			if (this->InfoOnCells)
			{
				// * lecture de ce qui reste, et on jete *
				UnusedResult=fscanf( stream, "%f", &tampon );
				UnusedResult=fscanf( stream, "%f", &tampon );
				UnusedResult=fscanf( stream, "%f", &tampon );
				UnusedResult=fscanf( stream, "%f", &tampon );
			}

		}
	
	// fermeture du fichier
	fclose( stream );
	
	}
	
	// copy de la liste de points dans la structure et desallocation de la memoire
	output_mesh->SetPoints(temp_points);
	temp_points->Delete();
	
	// copy de la liste de cellules dans la structure et desallocation de la memoire
	output_mesh->SetPolys(temp_cells);
	temp_cells->Delete();
}








