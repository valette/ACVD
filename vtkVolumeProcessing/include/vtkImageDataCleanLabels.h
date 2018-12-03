#include <vtkSimpleImageToImageFilter.h>



class VTK_EXPORT vtkImageDataCleanLabels : public vtkSimpleImageToImageFilter
{

public :

	static vtkImageDataCleanLabels *New();
	vtkTypeMacro(vtkImageDataCleanLabels,vtkSimpleImageToImageFilter);

protected :

	virtual void SimpleExecute (vtkImageData *input, vtkImageData *output);

	/// the constructor
	vtkImageDataCleanLabels()
	{
	}

	/// the destructor
	 ~vtkImageDataCleanLabels()
	{
	}
};

