/*=========================================================================

  Program:   vtkTag
  Module:    vtkSurface
  Language:  C++
  Date:      2008/04
  Auteur:    Sebastien VALETTE

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

#ifndef __vtkTag_h
#define __vtkTag_h

#include <limits.h>
#include "vtkObjectFactory.h"

class VTK_EXPORT vtkTag : public vtkObject
{

public:

	/// The Constructor vtkRecons::New();
	static vtkTag *New()
	{
		// First try to create the object from the vtkObjectFactory
		vtkObject *ret = vtkObjectFactory::CreateInstance ("vtkTag");
		if (ret)
		{
			return (vtkTag *) ret;
		}
		// If the factory was unable to create the object, then create it here.
		return (new vtkTag);
	}

	vtkTypeMacro(vtkTag,vtkObject);

	void SetNumberOfItems(int NumberOfItems)
	{
		this->LastTag->SetNumberOfValues(NumberOfItems);
		this->ResetCounter();
	}
	
	void Reset()
	{
		if (this->Time==INT_MAX)
			this->ResetCounter();
		else
			this->Time++;
	}
	
	bool IsTagged(vtkIdType Item)
	{
		return (this->LastTag->GetValue(Item)==this->Time);
	}
	
	void Tag(vtkIdType Item)
	{
		this->LastTag->SetValue(Item,this->Time);
	}

	void UnTag(vtkIdType Item)
	{
		this->LastTag->SetValue(Item,this->Time-1);
	}	

protected:

	void ResetCounter()
	{
		int NumberOfItems=this->LastTag->GetSize();
		this->Time=0;
		for (int i=0;i<NumberOfItems;i++)
			this->LastTag->SetValue(i,-1);
	}

	vtkIntArray *LastTag;
	
	int Time;

	/// constructor
	vtkTag()
	{
		this->LastTag=vtkIntArray::New();
	}

	/// desctructor
	virtual ~vtkTag()
	{
		this->LastTag->Delete();
	}
};

class VTK_EXPORT vtkTagWithList : public vtkTag
{

public:

	static vtkTagWithList *New()
	/// The Constructor vtkRecons::New();
	{
		// First try to create the object from the vtkObjectFactory
		vtkObject *ret = vtkObjectFactory::CreateInstance ("vtkTag");
		if (ret)
		{
			return (vtkTagWithList *) ret;
		}
		// If the factory was unable to create the object, then create it here.
		return (new vtkTagWithList);
	}

	vtkTypeMacro(vtkTagWithList,vtkObject);
	vtkGetObjectMacro(TaggedItems, vtkIdList)

	void Reset()
	{
		this->TaggedItems->Reset();
		this->vtkTag::Reset();
	}
	
	bool IsTagged(vtkIdType Item)
	{
		return (this->LastTag->GetValue(Item)==this->Time);
	}
	
	void Tag(vtkIdType Item)
	{
		if (!this->IsTagged(Item))
			this->TaggedItems->InsertNextId(Item);
		this->vtkTag::Tag(Item);
	}

protected:

	vtkIdList *TaggedItems;
	
	/// constructor
	vtkTagWithList()
	{
		this->TaggedItems=vtkIdList::New();
	}

	/// desctructor
	virtual ~vtkTagWithList()
	{
		this->TaggedItems->Delete();
	}
};

#endif
