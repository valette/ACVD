/*=========================================================================

Program:   sampleMesh : sample the surface of a mesh
Language:  C++
Date:      2015/10
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME sampleMesh
// .SECTION Description

#include <random>
#include <vtkTimerLog.h>
#include <vtkButterflySubdivisionFilter.h>
#include "vtkSurface.h"

int main( int argc, char *argv[] )
{
	if (argc < 3) {
		cout << "Usage : sampleMesh mesh numSamples [numberOfButterflySubdivisions]" << endl;
		exit(1);
	}

	vtkNew<vtkTimerLog> Timer;
	Timer->StartTimer();

	cout << "load : " << argv[ 1 ] << endl;
	vtkNew<vtkSurface> mesh;
	mesh->CreateFromFile( argv[ 1 ] );

	if ( argc > 3 ) {
		vtkNew<vtkButterflySubdivisionFilter> subdivision;
		subdivision->SetInputData( mesh );
		subdivision->SetNumberOfSubdivisions( atoi( argv[ 3 ] ) );
		subdivision->Update();
		mesh->CreateFromPolyData( subdivision->GetOutput() );
	}

	mesh->DisplayMeshProperties();
	int numTris = mesh->GetNumberOfCells();
	double sArea = 0;
	std::vector<float> sAreas;

	for (int i = 0; i < numTris; i++) {
		sArea += mesh->GetFaceArea(i);
		sAreas.push_back( sArea );
	}
	std::cout << "Total area : " << sArea << std::endl;

	int numSamples = atoi( argv[ 2 ] );
	vtkIdType v1, v2, v3;
	double P1[3], P2[3], P3[3], N[3], sample[3];
	std::random_device seeder;
	std::mt19937 generator( seeder() );
	std::uniform_real_distribution<float> rng( 0.0, sArea );
	std::uniform_real_distribution<float> rng2(0.0, 1.0);
	std::ofstream fileOutput;
	fileOutput.open ("output.xyz");

	for (int i = 0; i < numSamples; i++) {
		auto iter = lower_bound( sAreas.begin(), sAreas.end(), rng( generator ) );
		int index = iter - sAreas.begin();
		mesh->GetFaceVertices(index, v1, v2, v3);
		mesh->GetTriangleNormal(index, N);
		mesh->GetPointCoordinates(v1, P1);
		mesh->GetPointCoordinates(v2, P2);
		mesh->GetPointCoordinates(v3, P3);
		float r1 = rng2(generator);
		float r2 = rng2(generator);
		float x1 = 1 - sqrt( r1 );
		float x2 = sqrt( r1 ) * ( 1 - r2 );
		float x3 = r2 * sqrt( r1 );

		for (int k = 0; k < 3; k++)
			sample[k] = P1[k] * x1 + P2[k] * x2 + P3[k] * x3;

		fileOutput << sample[0] << '\t' << sample[1] << '\t'
			<< sample[2] << '\t' << N[0] << '\t' << N[1]
			<< '\t' << N[2];
		if ( i < numSamples - 1 ) fileOutput << std::endl;

	}
	Timer->StopTimer();
	fileOutput.close();
	cout << "Done in " << Timer->GetElapsedTime() << "s" << endl;
}
