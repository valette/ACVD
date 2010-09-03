/*=========================================================================

Module:	 vtkQuadricTools
Language:	 C++
Date:		 2008/12
Auteur:	Sebastien Valette,

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
#include "vtkQuadricTools.h"


// Adds to Quadric the quadric equivalent to a plane passing by the given Point with the Given Normal
// If FullQuadric=true then the whole quadric will be created (10 coefficients).
// Otherwise, only the 9 first coefficients will be given
void vtkQuadricTools::AddPointWithNormalQuadric(double *Quadric, double *Point,
	double *Normal, double Factor, bool FullQuadric)
{
	double d=-vtkMath::Dot(Normal,Point);
	Quadric[0]+=	Normal[0]*Normal[0]*Factor;
	Quadric[1]+=	Normal[0]*Normal[1]*Factor;
	Quadric[2]+=	Normal[0]*Normal[2]*Factor;
	Quadric[3]+=	Normal[0]*d*Factor;
	Quadric[4]+=	Normal[1]*Normal[1]*Factor;
	Quadric[5]+=	Normal[1]*Normal[2]*Factor;
	Quadric[6]+=	Normal[1]*d*Factor;
	Quadric[7]+=	Normal[2]*Normal[2]*Factor;
	Quadric[8]+=	Normal[2]*d*Factor;
	if (FullQuadric)
		Quadric[9]+=d*d*Factor;
}

// Adds the quadric computed from Face to Quadric, weighted by Factor
// If FullQuadric=true then the whole quadric will be created (10 coefficients).
// Otherwise, only the 9 first coefficients will be given
void vtkQuadricTools::AddTriangleQuadric(double *Quadric,vtkSurface *Mesh,int Face,bool FullQuadric)
{
	double x1[3],x2[3],x3[3];
	vtkIdType V1,V2,V3;
	double quadric4x4[4][4];

	Mesh->GetFaceVertices(Face,V1,V2,V3);
	Mesh->GetPointCoordinates(V1,x1);
	Mesh->GetPointCoordinates(V2,x2);
	Mesh->GetPointCoordinates(V3,x3);
	vtkTriangle::ComputeQuadric(x1,x2,x3,quadric4x4);

	Quadric[0]+=	quadric4x4[0][0];
	Quadric[1]+=	quadric4x4[0][1];
	Quadric[2]+=	quadric4x4[0][2];
	Quadric[3]+=	quadric4x4[0][3];
	Quadric[4]+=	quadric4x4[1][1];
	Quadric[5]+=	quadric4x4[1][2];
	Quadric[6]+=	quadric4x4[1][3];
	Quadric[7]+=	quadric4x4[2][2];
	Quadric[8]+=	quadric4x4[2][3];
	if (FullQuadric)
		Quadric[9]+=	quadric4x4[3][3];
}

// Computes the displacement needed to reach the best position according to the quadric
// MaxNumberOfUsedSingularValues defines the number of singular values used (generally 3)
int vtkQuadricTools::ComputeDisplacement(double *Quadric, double *Point, double *Displacement
			,int MaxNumberOfUsedSingularValues, double SVThreshold)
{
	int	i;
	double A[3][3], U[3][3], VT[3][3];
	double b[3], w[3];
	double tempMatrix[3][3];
	double tempMatrix2[3][3];

	int RankDeficiency=0;

	b[0] = -Quadric[3];
	b[1] = -Quadric[6];
	b[2] = -Quadric[8];

	A[0][0] = Quadric[0];
	A[0][1] = A[1][0]	= Quadric[1];
	A[0][2] = A[2][0]	= Quadric[2];
	A[1][1] = Quadric[4];
	A[1][2] = A[2][1]	= Quadric[5];
	A[2][2] = Quadric[7];

	vtkMath::SingularValueDecomposition3x3(A, U, w,	VT);

	// compute all eigen values absolute values
	double AbsolutesEigenValues[3];
	double maxW	= -1.0;
	for (int j=0;j<3;j++)
	{
		double AbsoluteEigenValue=fabs(w[j]);
		AbsolutesEigenValues[j]=AbsoluteEigenValue;
		if (AbsoluteEigenValue > maxW)
			maxW = AbsoluteEigenValue;
	}
	double invmaxW=1.0/maxW;

	for (i=0;i<3;i++)
	{
		double LocalMaxW=-1;
		int IndexMax=-1;

		// find the remaining eigenvalue with highest absolute value
		for (int j=0;j<3;j++)
		{
			if (LocalMaxW<AbsolutesEigenValues[j])
			{
				LocalMaxW=AbsolutesEigenValues[j];
				IndexMax=j;
			}
		}

		if (( AbsolutesEigenValues[IndexMax]*invmaxW > SVThreshold)
						&&(MaxNumberOfUsedSingularValues>0))
		{
			// If this is true,	then w[i] != 0,	so this	division is	ok.
			double Inv	=  1.0/w[IndexMax];
			tempMatrix[IndexMax][0]=U[0][IndexMax]*Inv;
			tempMatrix[IndexMax][1]=U[1][IndexMax]*Inv;
			tempMatrix[IndexMax][2]=U[2][IndexMax]*Inv;
		}
		else
		{
			tempMatrix[IndexMax][0]=0;
			tempMatrix[IndexMax][1]=0;
			tempMatrix[IndexMax][2]=0;
			RankDeficiency++;
		}

		// set the eigenvalu to -2 to remove it from subsequent tests
		AbsolutesEigenValues[IndexMax]=-2;
		MaxNumberOfUsedSingularValues--;
	}

	vtkMath::Transpose3x3(VT, VT);
	vtkMath::Multiply3x3(VT,	tempMatrix,	tempMatrix2);
	vtkMath::Multiply3x3(A,	Point,	Displacement);
	for	(i = 0;	i <	3; i++)
		Displacement[i] =	b[i] - Displacement[i];
	vtkMath::Multiply3x3(tempMatrix2, Displacement, Displacement);
	return (RankDeficiency);
}

// Projects the point on the position giving the minimum quadric error
// returns the rank deficiency of the quadric.
// MaxNumberOfUsedSingularValues defines the number of singular values used (generally 3)
int vtkQuadricTools::ComputeRepresentativePoint(double *Quadric, double *Point
				,int MaxNumberOfUsedSingularValues, double SVThreshold)
{
	double tempVector[3];
	int RankDefficiency=vtkQuadricTools::ComputeDisplacement(Quadric,Point,
							tempVector,MaxNumberOfUsedSingularValues,SVThreshold);
	for	(int i = 0;	i <	3; i++)
		Point[i] += tempVector[i];
	return RankDefficiency;
}

double vtkQuadricTools::Evaluate(double *Quadric, double *Point, bool FullQuadric)
{
	double Evaluation=Point[0]*Point[0]*Quadric[0]
	+Point[1]*Point[1]*Quadric[4]
	+Point[2]*Point[2]*Quadric[7]
	+2.0*Point[0]*Point[1]*Quadric[1]
	+2.0*Point[0]*Point[2]*Quadric[2]
	+2.0*Point[1]*Point[2]*Quadric[5]
	+2.0*Point[0]*Quadric[3]
	+2.0*Point[1]*Quadric[6]
	+2.0*Point[2]*Quadric[8];
	if (FullQuadric)
		Evaluation+=Quadric[9];
	return (Evaluation);
}

void vtkQuadricTools::GetPointQuadric(vtkSurface *Mesh, vtkIdType Vertex, double *Quadric, bool FullQuadric)
{
	for (int i=0;i<9;i++)
		Quadric[i]=0;
	
	if (FullQuadric)
		Quadric[9]=0;

	double Point[3];
	Mesh->GetPoint(Vertex,Point);

	Mesh->GetVertexNeighbourFaces(Vertex,this->List);
	for (vtkIdType i=0;i<this->List->GetNumberOfIds();i++)
		this->AddTriangleQuadric(Quadric, Mesh, this->List->GetId(i), FullQuadric);
}

vtkStandardNewMacro(vtkQuadricTools);

vtkQuadricTools::vtkQuadricTools()
{
	this->List=vtkIdList::New();
};

vtkQuadricTools::~vtkQuadricTools()
{
	this->List->Delete();
};	

