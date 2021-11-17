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

#include <vtkCommand.h>
#include "vtkSurface.h"
#include "RenderWindow.h"

void traverse( vtkSurface *test ) {

	vtkIdType v1,v2,v3;		// stockage de l'identite de points
	cout << "*******************************************" << endl;

	// display points coordinates
	for (int i=0;i<test->GetNumberOfPoints();i++)
	{
		double Point[3];
		test->GetPoint(i,Point);
		cout<<"Coordinates for vertex "<<i<<" : "<<Point[0]<<" "<<Point[1]<<" "<<Point[2]<<endl;
	}

	cout << "*******************************************" << endl;

	// loop on faces
	for (int i=0;i<test->GetNumberOfCells();i++)
	{
		test->GetFaceVertices(i,v1,v2,v3);
		cout<<"Face "<<i<<" has vertices "<<v1<<", "<<v2<<", "<<v3;
		if ( !test->IsFaceActive( i ) ) cout << " inactive";
		cout << endl;
	}

	cout << "*******************************************" << endl;

	// loop on edges
	for (int i=0;i<test->GetNumberOfEdges();i++)
	{
		test->GetEdgeVertices(i,v1,v2);
		cout<<"Edge "<<i<<" has vertices "<<v1<<" and "<<v2<<endl;
	}

	cout << "*******************************************" << endl;

	// loop on edges
	for (int i=0;i<test->GetNumberOfEdges();i++)
	{
		test->GetEdgeFaces(i,v1,v2);
		cout<<"Edge "<<i<<" has adjacent faces "<<v1<<" and "<<v2<<endl;
	}

	cout << "*******************************************" << endl;

}

// Here is a small sample program showing elementary methods to use vtkSurface

int main( int argc, char *argv[] )
{

	// points and cells definitions
	static double x[7][3]={{0,0,1},{0.707,0.707,0},{0.707,-0.707,0},{-0.707,-0.707,0},{-0.707,0.707,0},{0,0,-1}, {0,0,2}};
	static int pts[12][3]={{0,1,2},{0,2,3},{0,3,4},{0,4,1},{5,1,2},{5,2,3},{5,3,4},{5,4,1},{6,1,2},{6,2,3},{6,3,4},{6,4,1} }; 

	vtkSurface *test=vtkSurface::New();

	// create vertices
	for (int i=0;i<7;i++) test->AddVertex(x[i]);

	// create triangles
	for (int i=0;i<6;i++) 
		test->AddFace(pts[i][0],pts[i][1],pts[i][2]);

	traverse( test );

	// display mesh

	RenderWindow *Window=RenderWindow::New();
	Window->SetInputData(test);
	Window->Render();
	Window->Interact();

	// add more triangles
	for (int i=6;i<8;i++) 
		test->AddFace(pts[i][0],pts[i][1],pts[i][2]);

	traverse( test );

	//render again
	Window->Render();
	Window->Interact();

	// remove 4 first triangles
	for (int i=0;i<4;i++) test->DeleteFace( i );

	traverse( test );
	
	//render again
	Window->Render();
	Window->Interact();

	// add more triangles
	for (int i=8;i<12;i++) 
		test->AddFace(pts[i][0],pts[i][1],pts[i][2]);

	traverse( test );
	test->SQueeze();

	//render again
	Window->Render();
	Window->Interact();

	while( true ) {
		vtkIdType e = rand() % test->GetNumberOfEdges();
		cout << "flip edge : " << e << endl;
		test->FlipEdge(e);
		traverse( test );

		//render again
		Window->Render();
		Window->Interact();
	}

	//delete objects
	Window->Delete();
	test->Delete();

	return(0);
}
