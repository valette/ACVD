/*=========================================================================

  Program:   Mailleur 3D multi-résolution (Creatis 2000 ~ nowadays)
  Module:    vtkOFFWriter.cxx
  Language:  C++
  Date:      2003/05
  Auteurs:   Alexandre Gouaillard
=========================================================================*/

// .NAME vtkOFFWriter - write stereo lithography files
// .SECTION Description
// vtkOFFWriter writes OFF files (.off) files.

#ifndef __vtkOFFWriter_h
#define __vtkOFFWriter_h

#include <vtkPolyDataWriter.h>

class VTK_EXPORT vtkOFFWriter : public vtkPolyDataWriter
{
public:
  static vtkOFFWriter *New();
  vtkTypeMacro(vtkOFFWriter,vtkPolyDataWriter);
  //virtual void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkOFFWriter() {};
  ~vtkOFFWriter() {};

  void WriteData();

  private:
  vtkOFFWriter(const vtkOFFWriter&);  // Not implemented.
  void operator=(const vtkOFFWriter&);  // Not implemented.
};

#endif

