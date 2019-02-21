/*=========================================================================

Program:   sampleMesh : sample the surface of a mesh
Language:  C++
Date:      2015/10
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME sampleMesh
// .SECTION Description

#include <vtkTimerLog.h>
#include <vtkMinimalStandardRandomSequence.h>
#include "vtkSurface.h"


int main( int argc, char *argv[] )
{
	if (argc < 2) {
		cout << "Usage : sampleMesh mesh numSamples" << endl;
		exit(1);
	}

	vtkTimerLog *Timer = vtkTimerLog::New();
	Timer->StartTimer();

	cout << "load : " << argv[1] << endl;
	vtkSurface *mesh = vtkSurface::New();
	mesh->CreateFromFile(argv[1]);

	mesh->DisplayMeshProperties();
	int numTris = mesh->GetNumberOfCells();
	double sArea = 0;
	for (int i = 0; i < numTris; i++) {
		sArea += mesh->GetFaceArea(i);
	}
	std::cout << "Total area : " << sArea << std::endl;

	int numSamples = atoi(argv[2]);
	vtkIdType v1, v2, v3;
	double P1[3], P2[3], P3[3], N[3];
	int n2 = 0;
	double sample[3];

	vtkMinimalStandardRandomSequence *rng = vtkMinimalStandardRandomSequence::New();

	ofstream fileOutput;
	fileOutput.open ("output.xyz");

	for (int i = 0; i < numTris; i++) {
		double area = mesh->GetFaceArea(i);
		mesh->GetTriangleNormal(i,N);
		mesh->GetFaceVertices(i, v1, v2, v3);
		mesh->GetPointCoordinates(v1, P1);
		mesh->GetPointCoordinates(v2, P2);
		mesh->GetPointCoordinates(v3, P3);
		int n = floor(0.5 + numSamples * area / sArea);
		n2 +=n;
		for (int j = 0; j < n ; j++) {
			double x1 = rng->GetValue();
			rng->Next();
			double x2 = rng->GetValue();
			rng->Next();
			double x3 = rng->GetValue();
			rng->Next();
			double sum = x1 + x2 + x3;
			x1 /= sum;
			x2 /= sum;
			x3 /= sum;
			for (int k = 0; k < 3; k++) {
				sample[k] = P1[k] * x1 + P2[k] * x2 + P3[k] * x3;
			}
			fileOutput << sample[0] << '\t' << sample[1] << '\t'
				<< sample[2] << '\t' << N[0] << '\t' << N[1]
				<< '\t' << N[2] << std::endl;
		}
	}
	cout << n2 << " samples thrown" << endl;
	Timer->StopTimer();
	fileOutput.close();

	cout << "Done in " << Timer->GetElapsedTime() << "s" << endl;
}
