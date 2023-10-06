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
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkLight.h>
#include <vtkLightCollection.h>
#include <vtkPNGWriter.h>
#include <vtkBMPWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkPostScriptWriter.h>
#include <vtkSelectVisiblePoints.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkLabeledDataMapper.h>
#include <vtkActor2D.h>
#include <vtkCellArray.h>
#include <vtkIdFilter.h>
#include <vtkCellCenters.h>
#include <vtkTextProperty.h>
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkTextMapper.h>
#include <vtkCoordinate.h>
#include <vtkLight.h>
#include <vtkPolyDataNormals.h>
#include <vtkPLYWriter.h>
#include <vtkActor2DCollection.h>
#include <vtkSTLWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkOBJExporter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkRenderLargeImage.h>
#include <vtkTubeFilter.h>
#include <vtkTextProperty.h>
#include <vtkVRMLExporter.h>

#include "RenderWindow.h"
#include "vtkSurfaceIterators.h"

#define _RANDOM_SEED 1000

class MyInteractorStyleTrackballCamera:public
	vtkInteractorStyleTrackballCamera
{
public:
	RenderWindow * Window;

	static MyInteractorStyleTrackballCamera *New ()
	{
		// First try to create the object from the vtkObjectFactory
		vtkObject *ret =
			vtkObjectFactory::
			CreateInstance ("MyInteractorStyleTrackballCamera");
		if (ret)
		{
			return (MyInteractorStyleTrackballCamera *) ret;
		}
		// If the factory was unable to create the object, then create it here.
		return (new MyInteractorStyleTrackballCamera);
	}

	virtual void OnChar ()
	{
		vtkRenderWindowInteractor *rwi = this->Interactor;
		switch (rwi->GetKeyCode ())
		{
		case '1':
			this->Window->Capture ("Capture.jpg");
			return;
		case '2':
			this->Window->Capture ("Capture.eps");
			return;
		case '3':
			this->Window->Capture ("Capture.png");
			return;		
		case '4' :
			int Factor;
			cout<<"Enter new magnification factor :";
			cin>>Factor;
			this->Window->SetCaptureMagnificationFactor(Factor);
			return;	
		case 'E' :
			int Number;
			cout<<"Enter number of interactions to skip :";
			cin>>Number;
			this->Window->SkipInteractions(Number);
			return;
		case 'N':
		case 'n':
			this->Window->EnableNormalMap();
			return;
				
		case 'x':
		case 'X':
			if (this->EdgesDisplay == 0)
			{
				this->EdgesDisplay = 1;
				this->Window->DisplayInputEdges ();
			}
			else
			{
				this->EdgesDisplay = 0;
				this->Window->SetInputEdges (0);
			}
			this->Window->Render ();
			return;
		case 'I':
		case 'i':
			this->Window->SwitchDisplayIds ();
			return;
		case 'a':
			this->Window->Render();
			return;
		case 'O':
		case 'o':
			this->Window->SwitchOrientation();
			this->Window->Render();
			return;
		case '5':
		{
			char name[100];
			cout<<"Save Camera (.view)"<<endl;
			cout<<"Enter name  : ";
			cin>>name;

			this->Window->SaveCamera(name);

			return;
		}
		case '6':
		{
			char name[100];
			int flag;
			cout<<"Load Camera"<<endl;
			cout<<"Enter name  : ";
			cin>>name;

			flag = this->Window->LoadCamera(name);
			this->Window->Render();

			if(flag==1)
			cout<<"Load "<<name<<"	->	[OK]"<<endl;
			else
			cout<<"Load"<<name<<"	->	[BUG]"<<endl;
	
			return;
		}
		case '7':
		{
			vtkPolyDataWriter * SWriter = vtkPolyDataWriter::New ();
			SWriter->SetInputData (this->Window->GetInput ());
			SWriter->SetFileName ("mesh.vtk");
			SWriter->Write ();
			SWriter->Delete ();
			return;
		}		
		case '8':
		{
			vtkOBJExporter * SWriter = vtkOBJExporter::New ();
			SWriter->SetInput (this->Window->GetvtkRenderWindow());
			SWriter->SetFilePrefix ("mesh");
			SWriter->Write ();
			SWriter->Delete ();
			return;
		}
		case '9':
		{
			vtkSTLWriter * SWriter = vtkSTLWriter::New ();
			SWriter->SetInputData (this->Window->GetInput ());
			SWriter->SetFileName ("mesh.stl");
			SWriter->Write ();
			SWriter->Delete ();
			return;
		}
		case '0':
		{
			vtkPLYWriter * Writer = vtkPLYWriter::New ();
			Writer->SetInputData (this->Window->GetInput ());
			Writer->SetFileName ("mesh.ply");
			Writer->Write ();
			Writer->Delete ();
			return;
		}
		case '_':
		{
			vtkVRMLExporter * Writer = vtkVRMLExporter::New ();
			Writer->SetInput (this->Window->GetvtkRenderWindow());
			Writer->SetFileName ("mesh.wrl");
			Writer->Write ();
			Writer->Delete ();
			return;
		}
		case '+':
		{
			if (this->Window->GetEdgesActor())
			{
				this->Window->GetEdgesActor()->GetProperty()->SetLineWidth(
					this->Window->GetEdgesActor()->GetProperty()->GetLineWidth()+1);
				this->Window->Render();
			}
			return;
		}
		case '-':
		{
			if (this->Window->GetEdgesActor())
			{
				this->Window->GetEdgesActor()->GetProperty()->SetLineWidth(
					this->Window->GetEdgesActor()->GetProperty()->GetLineWidth()-1);
				this->Window->Render();
			}
			return;
		}
		case 'v':
		{
			vtkIdType Vertex;
			cout<<"Which Vertex do you want details from ? ";
			cin>>Vertex;
			double Point[3];
			this->Window->GetInput()->GetPoints()->GetPoint(Vertex,Point);
			cout<<"Coordinates : "<<Point[0]<<" "<<Point[1]<<" "<<Point[2]<<endl;
			/*
			vtkSurfaceVertexRingRandomIterator Iterator;
			Iterator.SetInput((vtkSurface *)this->Window->GetInput());
			cout<<"Which Vertex do you want details from ? ";
			cin>>Vertex;
			Iterator.InitTraversal(Vertex);
			cout<<"Neighbours :";
			while (1)
			{
				vtkIdType Neighbour=Iterator.GetNextVertex();
				if (Neighbour==-1)
					break;
				cout<<" "<<Neighbour;				
			}
			cout<<endl;
			Iterator.InitTraversal(Vertex);
			cout<<"Edges :";
			while (1)
			{
				vtkIdType Edge=Iterator.GetNextEdge();
				if (Edge==-1)
					break;
				cout<<" "<<Edge;				
			}
			cout<<endl;
			*/
			return;
		}
		case 't':
		{
			vtkIdType Face,v1,v2,v3;
			cout<<"Which triangle do you want details from? ";
			cin>>Face;
			((vtkSurface *)this->Window->GetInput())->GetFaceVertices(Face,v1,v2,v3);
			cout<<"Vertices : "<<v1<<" "<<v2<<" "<<v3<<endl;
		}
		default:
			cout<<"";
		}
		this->vtkInteractorStyleTrackballCamera::OnChar ();
	}

protected:
	MyInteractorStyleTrackballCamera ()
	{
		this->EdgesDisplay = 0;
	}
	~MyInteractorStyleTrackballCamera ()
	{
	};

private:
	int EdgesDisplay;

};

void RenderWindow::AttachToRenderWindow(RenderWindow *Window)
{
	if ( !Window || ( Window == this ) ) return;
	this->attachedWindows.insert( Window );
	std::set< RenderWindow* > allWindows;
	auto cam = Window->GetMeshRenderer()->GetActiveCamera();

	for ( auto win : this->attachedWindows )
		for ( auto attached : win->attachedWindows )
			allWindows.insert( attached );

	for ( auto win : this->attachedWindows ) {
		win->attachedWindows = allWindows;
		win->GetMeshRenderer()->SetActiveCamera( cam );
	}

}


void RenderWindow::Interact() 
{
	if (this->NumberOfInteractionsToSkip>0)
		this->NumberOfInteractionsToSkip--;
	else
		this->renWin->GetInteractor()->Start();
};

void RenderWindow :: SaveCamera(const char *filename)
{
	
	vtkCamera *camera=this->GetMeshRenderer()->GetActiveCamera();

	std::ofstream camera_file;
	camera_file.open (filename, std::ofstream::out | std::ofstream::trunc);
	double *ClippingRange = camera->GetClippingRange();
	double *FocalPoint = camera->GetFocalPoint();
	double *Position = camera->GetPosition();
	double *ViewUp = camera->GetViewUp();
	double ViewAngle = camera->GetViewAngle();
	camera_file<<"ClippingRange : "<<ClippingRange[0]<<" "<<ClippingRange[1]<<endl;
	camera_file<<"FocalPoint :"<<FocalPoint[0]<<" "<<FocalPoint[1]<<" "<<FocalPoint[2]<<endl;
	camera_file<<"Position : "<<Position[0]<<" "<<Position[1]<<" "<<Position[2]<<endl;
	camera_file<<"ViewUp : "<<ViewUp[0]<<" "<<ViewUp[1]<<" "<<ViewUp[2]<<endl;
	camera_file<<"ViewAngle : "<<ViewAngle<<endl;

	camera_file.close();
}

int RenderWindow :: LoadCamera(const char *filename)
{
	double ClippingRange[2];
	double FocalPoint[3];
	double Position[3];
	double ViewUp[3];
	double ViewAngle;

	vtkCamera *camera=this->GetMeshRenderer()->GetActiveCamera();
	std::ifstream file(filename);

   	if ( !file ) return 0;
	std::string line; // variable contenant chaque ligne lue
	std::string res;
	while ( std::getline( file, line) ) {
		//	std::cout<< line.c_str() << std::endl;
		sscanf(line.c_str(),"ClippingRange : %lf %lf",ClippingRange,ClippingRange+1);
		sscanf(line.c_str(),"FocalPoint : %lf %lf %lf",FocalPoint,FocalPoint+1,FocalPoint+2);
		sscanf(line.c_str(),"Position : %lf %lf %lf",Position,Position+1,Position+2);
		sscanf(line.c_str(),"ViewUp : %lf %lf %lf",ViewUp,ViewUp+1,ViewUp+2);
		sscanf(line.c_str(),"ViewAngle : %lf",&ViewAngle);
	}

	camera->SetClippingRange(ClippingRange);
	camera->SetFocalPoint(FocalPoint);
	camera->SetPosition(Position);
	camera->SetViewUp(ViewUp);
	camera->SetViewAngle(ViewAngle);
	file.close();
	
	return(1);
}


void RenderWindow::EnableNormalMap()
{
	vtkPolyDataNormals *Normals = vtkPolyDataNormals::New ();
	Normals->SetInputData (this->Input);
	Normals->SetFeatureAngle (60);
	Normals->Update();
	this->SetInputData (Normals->GetOutput ());
	Normals->Delete();
	this->Render();
}


void RenderWindow::SetCustomInteractorStyle()
{

	MyInteractorStyleTrackballCamera *CustomInteractorStyle=MyInteractorStyleTrackballCamera::New ();
	if (this->renWin->GetInteractor())
		this->renWin->GetInteractor()->SetInteractorStyle(CustomInteractorStyle);
	else
		cout<<"No Interactor"<<endl;
		
	CustomInteractorStyle->Window = this;
	CustomInteractorStyle->Delete();
}

void RenderWindow::SwitchOrientation()
{
	if (this->SInput) 
		this->SInput->SwitchOrientation();
}
	
void
RenderWindow::AddText (char *text, float x, float y, float z, float R,
		       float G, float B, int Size, int Justification)
{
	vtkTextMapper *Mapper;
	Mapper = vtkTextMapper::New ();

	vtkActor2D *Actor = vtkActor2D::New ();
	Mapper->SetInput (text);
	Mapper->GetTextProperty ()->SetColor (R, G, B);
	Mapper->GetTextProperty ()->SetFontSize (Size);
	Mapper->GetTextProperty ()->SetJustification(Justification);

	vtkCoordinate *coord = Actor->GetPositionCoordinate ();
	coord->SetCoordinateSystemToWorld ();
	coord->SetValue (x, y, z);

	Actor->SetMapper (Mapper);
	Mapper->Delete();
	this->GetMeshRenderer()->AddActor2D (Actor);
}

void
RenderWindow::SetText (const char *text)
{
	std::string s(text);
	if (s.length()==0)
	{
		// there is no text to display, remove the actor if present
		if (this->TextActor!=0)
		{
			this->GetMeshRenderer()->RemoveViewProp(this->TextActor);
// here the actor should be deleted, but this kills compiz on my machine... any hint?
//			this->TextActor->Delete();
			this->TextActor=0;
		}
	}
	else
	{
		if (this->TextActor==0)
		{
			vtkTextActor *Text=vtkTextActor::New();
			Text->SetTextScaleModeToViewport();
			Text->SetDisplayPosition(300, 50);
			Text->GetPosition2Coordinate()->SetCoordinateSystemToNormalizedViewport();
			this->TextActor=Text;
			vtkRenderer *Renderer=this->GetMeshRenderer();
			Renderer->AddActor(Text);
			vtkTextProperty *tprop = Text->GetTextProperty();
			tprop->SetFontSize(25);
			tprop->SetFontFamilyToArial();
			tprop->SetJustificationToCentered();
			tprop->BoldOn();
			tprop->SetColor(0, 0, 0);
		}
		this->TextActor->SetInput(text);
	}
}

vtkActor* RenderWindow::SetInputData (vtkPolyData * Input)
{
	vtkPolyDataMapper *Mapper=(vtkPolyDataMapper *) this->MeshActor->GetMapper();
	
	if ( Mapper==0 ) {
		Mapper=vtkPolyDataMapper::New();		
		Mapper->SetResolveCoincidentTopologyToPolygonOffset ();
		Mapper->SetRelativeCoincidentTopologyPolygonOffsetParameters(0,2);
		this->MeshActor->SetMapper(Mapper);
		Mapper->Delete();
	}

	Mapper->SetInputData(Input);
	if (this->lut) Mapper->SetLookupTable (this->lut);
	this->Input = Input;

//	this->MeshActor->GetProperty()->SetDiffuse(0.5);
//	this->MeshActor->GetProperty()->SetSpecular(0.5);
//	this->MeshActor->GetProperty()->SetColor(0.8,0.9,0);	

//	this->MeshActor->GetProperty()->SetAmbient(0.5);
//	this->MeshActor->GetProperty()->SetAmbientColor(0,0,0.9);
	return MeshActor;
}

vtkActor* RenderWindow::SetInputData (vtkSurface * Input)
{
	if (Input==0)
		return(0);
	this->SInput = Input;
	return (this->SetInputData((vtkPolyData *)Input));
}

void
RenderWindow::Capture (const char *filename)
{
	char fin[180];
	char FileName[180];
	const char *terminaison;
	double backuplinewidth=1 ;

	vtkRenderLargeImage *Capture =vtkRenderLargeImage::New ();
	Capture->SetInput (this->GetMeshRenderer());

	if (this->CaptureMagnificationFactor>1)
	{
		// we increase the edges width
		if (this->EdgesActor)
		{
					backuplinewidth=this->EdgesActor->GetProperty ()->GetLineWidth ();
					this->EdgesActor->GetProperty ()->SetLineWidth (backuplinewidth * CaptureMagnificationFactor);
		}
	}
	Capture->SetMagnification(this->CaptureMagnificationFactor);	

	strcpy (FileName, filename);
	if (FileName != NULL)
	{
		char *p;

		for (p = FileName; *p; ++p)
			*p = tolower (*p);
	}

	strcpy (fin, ".eps");
	terminaison = strstr (filename, fin);
	if (terminaison != NULL)
	{
		vtkPostScriptWriter *Writer = vtkPostScriptWriter::New ();
		Writer->SetInputConnection (Capture->GetOutputPort ());
		Writer->SetFileName (filename);
		Writer->Write ();
		Writer->Delete ();
	}
	strcpy (fin, ".bmp");
	terminaison = strstr (filename, fin);
	if (terminaison != NULL)
	{
		vtkBMPWriter *Writer = vtkBMPWriter::New ();
		Writer->SetInputConnection (Capture->GetOutputPort ());
		Writer->SetFileName (filename);
		Writer->Write ();
		Writer->Delete ();
	}
	strcpy (fin, ".png");
	terminaison = strstr (filename, fin);
	if (terminaison != NULL)
	{
		vtkPNGWriter *Writer = vtkPNGWriter::New ();
		Writer->SetInputConnection (Capture->GetOutputPort ());
		Writer->SetFileName (filename);
		Writer->Write ();
		Writer->Delete ();
	}
	strcpy (fin, ".jpg");
	terminaison = strstr (filename, fin);
	if (terminaison != NULL)
	{
		vtkJPEGWriter *Writer = vtkJPEGWriter::New ();
		Writer->SetInputConnection (Capture->GetOutputPort ());
		Writer->SetFileName (filename);
		Writer->SetQuality (85);
		Writer->Write ();
		Writer->Delete ();
	}

	if (this->CaptureMagnificationFactor>1)
	{
		// restore the edges width
		if (this->EdgesActor)
			this->EdgesActor->GetProperty ()->SetLineWidth (backuplinewidth);
	}

	Capture->Delete();
}

void
RenderWindow::SetLookupTable (vtkLookupTable * Colors)
{
	if ( this->lut && ( this->lut != Colors ) )
		this->lut->Delete();

	if ( Colors ) {

		this->lut = Colors;

	} else {

		this->lut=vtkLookupTable::New();
		double x[4];
		this->lut->SetRange (0.0, 2.0);
		this->lut->SetNumberOfTableValues (3);
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
		x[3] = 1;
		this->lut->SetTableValue (0, x);
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
		x[3] = 1;
		this->lut->SetTableValue (1, x);
		x[0] = 0;
		x[1] = 0;	
		x[2] = 1;
		x[3] = 0.1;
		this->lut->SetTableValue (2, x);	
	}
	
	if ( this->Input ) {

		vtkPolyDataMapper *Mapper=(vtkPolyDataMapper *) this->MeshActor->GetMapper();
		Mapper->SetLookupTable (this->lut );
		Mapper->SetScalarRange (this->lut ->GetTableRange ());

	}

}

void
RenderWindow::DisplayVerticesColors (vtkIntArray * Scalars)
{
	if (Scalars)
		this->Input->GetPointData ()->SetScalars (Scalars);

	int i;
	int max = 0;

	for (i = 0; i < this->Input->GetNumberOfPoints (); i++)
	{
		if (max < Scalars->GetValue (i))
			max = Scalars->GetValue (i);
	}
	this->DisplayRandomColors(max);	
}

void
RenderWindow::DisplayCellsColors (vtkIntArray * Scalars)
{
	if (Scalars)
		this->Input->GetCellData ()->SetScalars (Scalars);

	int i;
	int max = 0;

	for (i = 0; i < this->Input->GetNumberOfCells (); i++)
	{
		if (max < Scalars->GetValue (i))
			max = Scalars->GetValue (i);
	}
	this->DisplayRandomColors(max);	
}

void
RenderWindow::DisplayRandomColors (int NumberOfColors)
{
	vtkMath *Math = vtkMath::New ();
	Math->RandomSeed (_RANDOM_SEED);
	double x[4];

	vtkLookupTable *lut=this->lut;
	if ( !lut ) {
		lut = vtkLookupTable::New();
		lut->SetNumberOfTableValues( 0 );
	}
	lut->SetRange (0.0, NumberOfColors - 1);
	int oldSize = lut->GetNumberOfTableValues();
	if ( oldSize > 0 ) oldSize--;
	lut->SetNumberOfTableValues (NumberOfColors);

	for (int i = oldSize; i < NumberOfColors - 1; i++)
	{
		x[0] = x[1] = x[2] = Math->Random ();
		x[1] = Math->Random ();
		x[2] = Math->Random ();
		x[3] = 1;
		lut->SetTableValue (i, x);
	}

	x[0] = 1.0;
	x[1] = 1.0;
	x[2] = 1.0;
	x[3] = 1;
	lut->SetTableValue (NumberOfColors - 1, x);
	this->SetLookupTable(lut);
	Math->Delete ();
}

void RenderWindow::HighLightEdges(vtkIdList *Edges, double Radius)
{
	if(Edges!=0)
	{
		vtkPolyData *EdgesP=vtkPolyData::New();
		vtkIdList *vl=vtkIdList::New();
		vl->SetNumberOfIds(2);

		EdgesP->SetPoints(this->SInput->GetPoints());
		EdgesP->Allocate(Edges->GetNumberOfIds());
		for (int i=0;i<Edges->GetNumberOfIds();i++)
		{
			vtkIdType v1,v2;
			this->SInput->GetEdgeVertices(Edges->GetId(i),v1,v2);
			vl->SetId(0,v1);
			vl->SetId(1,v2);
			EdgesP->InsertNextCell(VTK_LINE,vl);
		}

		EdgesP->Modified();
		vl->Delete();
		
		vtkTubeFilter *Tube=vtkTubeFilter::New();
		Tube->SetRadius(Radius);
		Tube->SetNumberOfSides(20);
		Tube->SetInputData(EdgesP);

		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New( );
		mapper->SetInputData( Tube->GetOutput( ) );

		if( !this->HighlightedEdgesActor )
			HighlightedEdgesActor = vtkActor::New( );
		else
			this->GetMeshRenderer()->RemoveActor( HighlightedEdgesActor );

		HighlightedEdgesActor->SetMapper( mapper );
		HighlightedEdgesActor->GetProperty( )->SetColor( 1., 0., 0. );

		this->GetMeshRenderer()->AddActor( HighlightedEdgesActor );
		mapper->Delete( );
		Tube->Delete();
		EdgesP->Delete();
	}
	else 
		if( this->HighlightedEdgesActor )
			this->GetMeshRenderer()->RemoveActor( HighlightedEdgesActor );
}

void RenderWindow::HighLightVertices(vtkIdList *Vertices, double Radius)
{
	if(Vertices!=0)
	{
		double p[3];
		vtkPoints* points = vtkPoints::New( );

		for( vtkIdType i = 0; i < Vertices->GetNumberOfIds(); i++ )
		{
			this->SInput->GetPointCoordinates( Vertices->GetId(i), p );
			points->InsertNextPoint( p );            
		}

		vtkUnstructuredGrid* ps = vtkUnstructuredGrid::New( );
		ps->SetPoints( points );

		vtkSphereSource* s_sphere = vtkSphereSource::New( );
		s_sphere->SetRadius( Radius );
		int Res=20;
		s_sphere->SetThetaResolution(Res);
		s_sphere->SetPhiResolution(Res);
		s_sphere->Update( );

		vtkGlyph3D* glyph = vtkGlyph3D::New( );
		glyph->SetInputData( ps );
		glyph->SetSourceData( s_sphere->GetOutput( ) );
		glyph->Update( );

		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New( );
		mapper->SetInputData( glyph->GetOutput( ) );

		if( !this->HighlightedVerticesActor )
			HighlightedVerticesActor = vtkActor::New( );
		else
			this->GetMeshRenderer()->RemoveActor( HighlightedVerticesActor );

		HighlightedVerticesActor->SetMapper( mapper );
		HighlightedVerticesActor->GetProperty( )->SetColor( 1., 0., 0. );

		this->GetMeshRenderer()->AddActor( HighlightedVerticesActor );

		ps->Delete( );
		s_sphere->Delete( );
		glyph->Delete( );
		mapper->Delete( );
		points->Delete( );
	}
	else 
		if( this->HighlightedVerticesActor )
			this->GetMeshRenderer()->RemoveActor( HighlightedVerticesActor );
}

void RenderWindow::DisplayNonManifoldVertices( double Radius )
{
	if (this->SInput==0)
	{
		cout<<"Warning : trying to display non manifold vertices does not work for vtkPolyData inputs"<<endl;
		cout<<"Use vtkSurface to use this feature"<<endl;
		return;
	}

	vtkIdList *Vertices=vtkIdList::New();
    for( vtkIdType i = 0; i < this->SInput->GetNumberOfPoints( ); i++ )
    {
        if( !this->SInput->IsVertexManifold( i ) )
        	Vertices->InsertNextId(i);
    }
    this->HighLightVertices(Vertices,Radius);
    Vertices->Delete();
}

vtkActor *RenderWindow::AddPolyData (vtkPolyData * Input)
{

	vtkPolyDataMapper *mapper = vtkPolyDataMapper::New ();
	mapper->SetResolveCoincidentTopologyToPolygonOffset ();
	mapper->SetRelativeCoincidentTopologyPolygonOffsetParameters(0,2);
//	mapper->SetResolveCoincidentTopologyToShiftZBuffer();

	mapper->SetInputData (Input);

	vtkActor *Mactor = vtkActor::New ();
	Mactor->SetMapper (mapper);
	if (this->lut)
	{
		mapper->SetLookupTable (this->lut);
		mapper->SetScalarRange (this->lut->GetTableRange ());
	}

	this->GetMeshRenderer()->AddActor (Mactor);
//	Mactor->GetProperty ()->EdgeVisibilityOn ();
	Mactor->GetProperty ()->SetPointSize (7);
	Mactor->GetProperty()->SetEdgeColor  (0,0,0);


	Mactor->Delete ();
	mapper->Delete ();
	return (Mactor);
}

void
RenderWindow::DisplayInputEdges ()
{
	if (!this->SInput)
	{
		if (!this->Input)
			return;
			
		// The input is a PolyData : we have to manually build the edge list.
		int i,j,NPoints=this->Input->GetNumberOfPoints();
		vtkEdgeTable *Table=vtkEdgeTable::New();
		Table->InitEdgeInsertion(NPoints);
		vtkIdList *VList=vtkIdList::New();
		this->Input->BuildCells();
		for (i=0;i<this->Input->GetNumberOfCells();i++)
		{
			this->Input->GetCellPoints(i,VList);
			for (j=0;j<VList->GetNumberOfIds();j++)
			{
				if (Table->IsEdge(VList->GetId(j), VList->GetId((j+1)%VList->GetNumberOfIds()))<0)
					Table->InsertEdge(VList->GetId(j), VList->GetId((j+1)%VList->GetNumberOfIds()));
			}			
		}
		
		vtkIdType v1,v2;
		VList->SetNumberOfIds(2);

		vtkPolyData *Edges=vtkPolyData::New();

		Edges->Allocate(Table->GetNumberOfEdges());
		Edges->SetPoints(this->Input->GetPoints());
		Table->InitTraversal();
		for (i=0;i<Table->GetNumberOfEdges();i++)
		{
			Table->GetNextEdge(v1,v2);
			VList->SetId(0,v1);
			VList->SetId(1,v2);
			Edges->InsertNextCell(VTK_LINE,VList);
		}
		this->SetInputEdges (Edges);
		Edges->Delete();
		VList->Delete();
		Table->Delete();
		return;
		
	}

	vtkPolyData *Edges = this->SInput->GetEdgesPolyData ();
	this->SetInputEdges (Edges);
	Edges->Delete ();
}

void
RenderWindow::SetInputEdges (vtkPolyData * Edges)
{
	if (this->EdgesActor)
	{
			this->GetMeshRenderer()->RemoveActor (EdgesActor);
			this->EdgesActor=0;
	}
	
	if (Edges)
	{
	
		vtkPolyDataMapper *EdgesMapper = vtkPolyDataMapper::New ();
		EdgesMapper->SetResolveCoincidentTopologyToPolygonOffset ();
		EdgesMapper->SetRelativeCoincidentTopologyPolygonOffsetParameters(0,2);

		this->EdgesActor = vtkActor::New ();
		vtkIntArray *EdgesColor = vtkIntArray::New ();

		EdgesMapper->SetInputData (Edges);
		int i;

		EdgesColor->SetNumberOfValues (Edges->GetNumberOfCells ());
						     
		for (i = 0; i < Edges->GetNumberOfCells (); i++)
			EdgesColor->SetValue (i, 1);

		Edges->GetCellData ()->SetScalars (EdgesColor);
		EdgesColor->Delete();
		this->EdgesActor->GetProperty ()->SetLineWidth (1);
		this->EdgesActor->SetMapper (EdgesMapper);
		EdgesMapper->Delete();

		double x[4];
		vtkLookupTable *lut=vtkLookupTable::New();
		lut->SetRange (0.0, 3.0);
		lut->SetNumberOfTableValues (4);
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
		x[3] = 0;
		lut->SetTableValue (0, x);
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
		x[3] = 0.5;
		lut->SetTableValue (1, x);
		x[0] = 0;
		x[1] = 0;
		x[2] = 1;
		x[3] = 0.1;
		lut->SetTableValue (2, x);
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
		x[3] = 1;
		lut->SetTableValue (3, x);		
		EdgesMapper->SetLookupTable (lut);
		EdgesMapper->SetScalarRange (lut ->GetTableRange ());
		lut->Delete();
		this->GetMeshRenderer()->AddActor (EdgesActor);
		EdgesActor->Delete();
		
		EdgesActor->SetOrientation(this->MeshActor->GetOrientation());

	}
}

vtkIntArray*
RenderWindow::GetEdgesVisibilityArray()
{
	if (this->EdgesActor)
	{
		vtkPolyDataMapper *Mapper=(vtkPolyDataMapper *) this->EdgesActor->GetMapper();
		return ((vtkIntArray *) Mapper->GetInput()->GetCellData()->GetScalars());
	}
	else
		return 0;
}

void
RenderWindow::SwitchDisplayIds ()
{
	vtkActor2DCollection *Actors = this->GetMeshRenderer()->GetActors2D ();

	vtkActor2D *Actor;
	Actors->InitTraversal ();
	Actor = Actors->GetNextActor2D ();

	while (Actor)
	{
		if (Actor == this->rectActor)
		{
			// Display was On.
			this->SetDisplayIdsOff ();
			return;
		}
		Actor = Actors->GetNextActor2D ();

	}

	// display was Off
	this->SetDisplayIdsOn ();
}

void
RenderWindow::SetDisplayIdsOff ()
{

	vtkRenderer *Renderer=this->GetMeshRenderer();
	
	if (!this->rectActor)
		return;

	Renderer->RemoveActor (rectActor);
	Renderer->RemoveActor (pointLabels);
	Renderer->RemoveActor (cellLabels);
	
	rectActor->Delete();
	rectActor=0;
	pointLabels->Delete();
	cellLabels->Delete();

	this->Render ();
}


void
RenderWindow::SetDisplayIdsOn ()
{
	vtkRenderer *Renderer=this->GetMeshRenderer();

	if (!this->rectActor)
	{

		// create the actors
		this->rectActor = vtkActor2D::New ();
		this->pointLabels = vtkActor2D::New ();
		this->cellLabels = vtkActor2D::New ();


		//Begin of Visualize Point Ids and Cell Ids

		int xmin = 0;
		int ymin = 0;
		int xmax=this->renWin->GetSize()[0];
		int ymax=this->renWin->GetSize()[1];

		//Create sub-viewport window
		vtkPoints *pts = vtkPoints::New ();
		pts->InsertPoint (0, xmin, ymin, 0);
		pts->InsertPoint (1, xmax, ymin, 0);
		pts->InsertPoint (2, xmax, ymax, 0);
		pts->InsertPoint (3, xmin, ymax, 0);

		vtkCellArray *rect = vtkCellArray::New ();
		rect->InsertNextCell (5);
		rect->InsertCellPoint (0);
		rect->InsertCellPoint (1);
		rect->InsertCellPoint (2);
		rect->InsertCellPoint (3);
		rect->InsertCellPoint (0);

		vtkPolyData *selectRect = vtkPolyData::New ();
		selectRect->SetPoints (pts);
		selectRect->SetLines (rect);

		vtkPolyDataMapper2D *rectMapper = vtkPolyDataMapper2D::New ();
		rectMapper->SetInputData (selectRect);

		rectActor->SetMapper (rectMapper);
		vtkDataSet *DataSet = this->GetMeshActor()->GetMapper()->GetInput();

		//Generate ids for labeling
		vtkIdFilter *ids = vtkIdFilter::New ();
		ids->SetInputData (DataSet);
		ids->PointIdsOn ();
		ids->CellIdsOn ();
		ids->FieldDataOn ();
		ids->Update ();

		//Create labels for points
		vtkSelectVisiblePoints *visPts =
			vtkSelectVisiblePoints::New ();
		visPts->SetInputData (ids->GetOutput ());
		visPts->SetRenderer (Renderer);
		visPts->SelectionWindowOn ();

		visPts->SetSelection (xmin, xmax, ymin, ymax);
		visPts->SelectInvisibleOff();
		vtkLabeledDataMapper *ldm = vtkLabeledDataMapper::New ();
		ldm->SetInputConnection ( visPts->GetOutputPort ());
//		ldm->SetLabelModeToLabelIds ();
		ldm->SetLabelModeToLabelFieldData ();

		pointLabels->SetMapper (ldm);

		//Create labels for cells
		vtkCellCenters *cc = vtkCellCenters::New ();
		cc->SetInputData (ids->GetOutput ());

		vtkSelectVisiblePoints *visCells =
			vtkSelectVisiblePoints::New ();
		visCells->SetInputConnection (cc->GetOutputPort ());
		visCells->SetRenderer (this->GetMeshRenderer());
		visCells->SelectInvisibleOff();
		visCells->SelectionWindowOn ();
		visCells->SetSelection (xmin, xmax, ymin, ymax);

		vtkLabeledDataMapper *cellMapper =
			vtkLabeledDataMapper::New ();
		cellMapper->SetInputConnection (visCells->GetOutputPort ());
//		cellMapper->SetLabelFormat ("%g");
		cellMapper->SetLabelModeToLabelFieldData ();
//		cellMapper->SetLabelModeToLabelIds ();
		cellMapper->GetLabelTextProperty ()->SetColor (0, 1, 0);


		cellLabels->SetMapper (cellMapper);
	}

	Renderer->AddActor2D (rectActor);
	Renderer->AddActor2D (pointLabels);
	Renderer->AddActor2D (cellLabels);
	

	renWin->Render ();

}

RenderWindow *
RenderWindow::New ()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject *ret = vtkObjectFactory::CreateInstance ("RenderWindow");
	if (ret)
	{
		return (RenderWindow *) ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new RenderWindow);

}

RenderWindow::RenderWindow () : HighlightedVerticesActor( 0 ),HighlightedEdgesActor (0)	//constructeur
{
	vtkRenderer *Renderer=vtkRenderer::New();
	
	this->MeshActor = vtkActor::New ();
	this->MeshActor->GetProperty ()->SetEdgeColor (0, 0, 0);
	this->MeshActor->GetProperty ()->SetPointSize (10);
	
//     this->Renderer->SetBackground( 0.1, 0.2, 0.4 );
	Renderer->SetBackground (1.0, 1.0, 1.0);
	vtkLight *Light = vtkLight::New ();
	Light->SetLightTypeToCameraLight ();
	Light->SetPosition (0.2, 0.2, 1);
	
	Renderer->AddLight (Light);
	Light->Delete();

	this->renWin = vtkRenderWindow::New ();

	this->renWin->AddRenderer (Renderer);
	Renderer->Delete();

//      this->renWin->SetSize( 320, 700 ); // happy budha
	this->renWin->SetSize (600, 600);
//      this->renWin->SetSize( 600, 400 ); // dragon

	this->renWin->SetInteractor(vtkRenderWindowInteractor::New());
//	this->renWin->GetInteractor()->Delete();
	this->SetCustomInteractorStyle();
	vtkNew<vtkCallbackCommand> getOrientation;
	getOrientation->SetCallback(InteractionCallback);
	getOrientation->SetClientData(this);
	this->renWin->GetInteractor()->AddObserver(vtkCommand::InteractionEvent, getOrientation);
	this->attachedWindows.insert( this );

	this->lut = 0;

	this->Input = 0;
	this->SInput = 0;

	this->EdgesActor = 0;

	this->CaptureMagnificationFactor = 1;

	Renderer->AddActor (MeshActor);
	MeshActor->Delete();

	this->rectActor = 0;
	this->pointLabels = 0;
	this->cellLabels = 0;
	this->NumberOfInteractionsToSkip=0;
	this->TextActor=0;
}

RenderWindow::~RenderWindow ()	//Destructeur
{
	if (this->renWin)
		this->renWin->Delete ();

	if (this->lut)
		this->lut->Delete ();

	if (this->rectActor)
		this->rectActor->Delete ();

	if (this->pointLabels)
		this->pointLabels->Delete ();

	if (this->cellLabels)
		this->cellLabels->Delete ();
}
