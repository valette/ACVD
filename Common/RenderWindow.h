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

#ifndef _RENDERWINDOW_H_
#define _RENDERWINDOW_H_

#include <set>

#include <vtkDataSetMapper.h>
#include <vtkRenderWindow.h>
#include <vtkCommand.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkTextActor.h>

#include "vtkLookupTable.h"
#include "vtkSurface.h"

/**
 * An efficient class for 3D triangular mesh display.
 * It can display vtkPolyData and vtkSurface
 * 
 *  Custom interactions:
 *  - "0"			Save input mesh to "mesh.ply"
 *  - "9"			Save input mesh to "mesh.stl"
 *  - "8"			Save input mesh to "mesh.obj"
 *  - "7"			Save input mesh to "mesh.vtk"
 *  - "6"			Load camera viewport from file
 *  - "5"			Save camera viewport to file
 *
 *  - "1"			Capture window to "capure.jpg" image file
 *  - "2"			Capture window to "capture.eps" image file
 *  - "3"			Capture window to "capture.png" image file
 *  - "4"			Set the window size multiplication factor (>=1) to create high quality capture images	
 *
 *  - "E"			Set the number of times Interact() will be skipped (usefull for debuging)
 *  - "x"			Toggles On/off the display of the mesh edges on the surface
 *  - "n"			Enable Normal Mapping (displays smooth meshes)
 *  - "i"			Toggles On/Off the display of vertices and cells Ids in a subwindow
 *  - "a"			Forces Rendering
 *  - "+"			Increases the width of rendered edges
 *  - "-"			Decreases the width of rendered edges
 */

class  VTK_EXPORT RenderWindow
{
public:

	// The constructor, vtk style (although RenderWindow is not a vtkObject anymore...)
	static RenderWindow *New();
	void Delete() {delete this;};

	// Adds some text in the window,
	// The x,y and z coordinates are world coordinates, so the text will move with the object
	// R, G and B are the colors components
	// Size is the text Size
	void AddText(char *text,float x, float y, float z, float R, float G, float B, int Size, int Justification=1);
	
	// write 2D text on the screen;
	void SetText (const char *text);

	/// captures the Window to a file
	/// supported file formats are: BMP, EPS, PNG (automatic file format detection)
	void Capture(const char *filename);

	/// Sets the magnification factor for the capture (the dimensions of the image will be multiplied by Factor)
	void SetCaptureMagnificationFactor (int Factor){this->CaptureMagnificationFactor=Factor;};

	static void InteractionCallback(vtkObject* caller, long unsigned int evId,
						 void* clientData, void* /*callData*/)
	{
		RenderWindow *self = reinterpret_cast<RenderWindow*>(clientData);
		for ( auto win : self->attachedWindows ) {
			if ( win != self ) win->Render();
		}
	};

	/// Links the window viewport to an other RenderWindow
	void AttachToRenderWindow(RenderWindow *Window);

	/// returns the vtkCamera of the Window (usefull to change the view : rotations, zoom etc)
	vtkCamera *GetCamera() {return this->GetMeshRenderer()->GetActiveCamera();}

	/// Saves the current camera view to file (only the camera, not the object)
	void SaveCamera(const char *filename);
	
	/// Loads the camera view from file
	int LoadCamera(const char *filename);

	/// Returns the Renderer used in the window
	vtkRenderer *GetMeshRenderer()
	{
		vtkRendererCollection *Renderers=this->renWin->GetRenderers();
		Renderers->InitTraversal();
		return (Renderers->GetNextItem());
	}

	/// sets the size of the window in pixels
	void SetSize (int x, int y) {this->renWin->SetSize( x,y );};

	/// renders the scene
	void Render() {this->renWin->Render();};

	/// to add a PolyData to the scene
	vtkActor* AddPolyData(vtkPolyData *Input);

	/// Starts interactive rendering (viewport control with the mouse)
	void Interact();
	
	/// Set the number of times a call to RenderWindow::Interact() will be skipped
	void SkipInteractions (int NumberOfInteractionsToSkip)
	{ this->NumberOfInteractionsToSkip=NumberOfInteractionsToSkip;};

	/// Sets the name of the window
	void SetWindowName(const char *Name){this->renWin->SetWindowName(Name);};
	
	/// Sets the input PolyData
	vtkActor* SetInputData(vtkPolyData *Input);

	/// Sets the input VtkSurface
	vtkActor* SetInputData(vtkSurface *Input);

	/// Returns the input mesh (Warning: its is returned as a vtkPolyData but it may be a vtkSurface!)
	vtkPolyData *GetInput()
	{
		if (this->SInput)
			return (this->SInput);
		else
			return (this->Input);
	}

	/// Displays the Ids of the vertices and cells into a sub-window
	void SetDisplayIdsOn();
	
	/// stops displaying the Ids of the vertices and cells
	void SetDisplayIdsOff();
	
	/// switch between Off/On modes for Ids display
	void SwitchDisplayIds();

	/// Displays the mesh with the edges. Does not keep track of connectivity modification,
	/// so whenever the input mesh connectivity is modified, this method has to be called again
	void DisplayInputEdges();

	/// Adds a set of custom edges to the window
	void SetInputEdges(vtkPolyData *Edges);
	
	// Enables Normal Mapping (to display smooth meshes)
	void EnableNormalMap();
	
	/// Returns the vtkIntArray defining which edges are visible
	vtkIntArray* GetEdgesVisibilityArray();

	/// replace the default LookUpTable by a custom one
	void SetLookupTable(vtkLookupTable *Colors=0);
	
	/// returns the used vtkLoopupTable for coloring
	vtkLookupTable *GetLookupTable() {return (this->lut);};

	/// Displays the colors of the vertices given the Scalars array, with random coloring
	void DisplayVerticesColors(vtkIntArray *Scalars);

	/// Displays the colors of the cells given the Scalars array, with random coloring
	void DisplayCellsColors(vtkIntArray *Scalars);

	/// Creates random colors for rendering
	void DisplayRandomColors(int NumberOfColors);

	/// Puts spheres on each Vertex present on the list (Radius is the radius of the sphere)
	void HighLightVertices(vtkIdList *Vertices, double Radius);

	/// Puts tubes on each edge present on the list (Radius is the radius of the tube)
	void HighLightEdges(vtkIdList *Edges, double Radius);

	/// Puts spheres on each non-manifold vertex (Radius is the radius of the sphere)
    void DisplayNonManifoldVertices(double Radius);

	/// returns the Actor containing the mesh edges
	vtkActor *GetEdgesActor() {return (this->EdgesActor);};
	
	vtkActor *GetMeshActor() {return (this->MeshActor);};

	/// returns the underlying vtkRenderWindow
	vtkRenderWindow *GetvtkRenderWindow() {return (this->renWin);};
	
	/// Switches the orientation of the input mesh	
	void SwitchOrientation();
	
	// method to set the interactor style to allow additionnal function keys
	virtual void SetCustomInteractorStyle();

protected:

	// The window
	vtkRenderWindow *renWin;

	// The parameter defining how many times Interaction will be skiped (usually 0)	
	int NumberOfInteractionsToSkip;

	// The actor for the mesh
	vtkActor *MeshActor;
	
	// The Actor for the input edges
	vtkActor *EdgesActor;

    vtkActor *HighlightedVerticesActor;
    vtkActor *HighlightedEdgesActor;
	
	// The input mesh (when it is a vtkPolyData)
	vtkPolyData *Input;

	// The input mesh (when it is a vtkSurface)
	vtkSurface *SInput;
	
	// The lookup table used to store colors
	vtkLookupTable *lut;

	// the window dimensions will be multiplied by this value when saving a snapshot
	int CaptureMagnificationFactor;

	// 2D actors used to display the IDs on screen.
	vtkActor2D *rectActor;
	vtkActor2D *pointLabels;
	vtkActor2D *cellLabels;

	vtkTextActor *TextActor;

	std::set< RenderWindow* > attachedWindows;

	RenderWindow(); 
	virtual ~RenderWindow();

};

#endif
