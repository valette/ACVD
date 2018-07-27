/*=========================================================================

  Program:   Visualization Toolkit
  Language:  C++
  Thanks:    Thanks to Abdalmajeid M. Alyassin who developed this class.



=========================================================================*/
// .NAME vtkVolumeProperties - estimate volume, area, shape index of triangle mesh
// .SECTION Description
// vtkVolumeProperties estimates the volume, the surface area, and the
// normalized shape index of a model.  The algorithm implemented here is
// based on the discrete form of the divergence theorem.  The general
// assumption here is that the model is of closed surface.  For more
// details see the following reference (Alyassin A.M. et al, "Evaluation
// of new algorithms for the interactive measurement of surface area and
// volume", Med Phys 21(6) 1994.).
// NOTE: currently only triangles are processed. Use vtkTriangleFilter
// to convert any strips or polygons to triangles.


#ifndef __vtkVolumeProperties_h
#define __vtkVolumeProperties_h

#include <vtkPolyData.h>

class VTK_EXPORT vtkVolumeProperties : public vtkObject
{
public:
  // Description:
  // Constructs with initial values of zero.
  static vtkVolumeProperties *New();

  vtkTypeMacro(vtkVolumeProperties,vtkObject);
  //void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Compute and return the volume.
  double GetVolume() { this->Update();return this->Volume;}
  double GetSignedVolume() { this->Update();return this->SignedVolume;}

  // Description:
  // Compute and return the barycentric coordinates
  double GetXG() { this->Update();return this->XG;}
  double GetYG() { this->Update();return this->YG;}
  double GetZG() { this->Update();return this->ZG;}

  // Description:
  // Compute and return the area.
  double GetSurfaceArea() { this->Update();return this->SurfaceArea;}

  // Description:
  // Compute and return the normalized shape index. This characterizes the
  // deviation of the shape of an object from a sphere. A sphere's NSI
  // is one. This number is always >= 1.0.
  double GetNormalizedShapeIndex() { this->Update();return this->NormalizedShapeIndex;}

  void Execute();
  void Update();
  
  vtkSetObjectMacro(InputData, vtkPolyData);
  vtkGetObjectMacro(InputData, vtkPolyData);

protected:
  vtkVolumeProperties();
  ~vtkVolumeProperties();
  vtkVolumeProperties(const vtkVolumeProperties&) {};
  void operator=(const vtkVolumeProperties&) {};

  vtkPolyData *InputData;

  double  SurfaceArea;
  double  Volume;
  double  SignedVolume;
  double  XG;
  double  YG;
  double  ZG;
  double  NormalizedShapeIndex;
  vtkTimeStamp ExecuteTime;

};

#endif


