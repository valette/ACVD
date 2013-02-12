#include <vector>
#include <map>
#include <queue>
#include <vtkObjectFactory.h>
#include <vtkImageData.h>
#include "vtkImageDataCleanLabels.h"

vtkStandardNewMacro(vtkImageDataCleanLabels);


struct voxel {
	int pointer;
	int x;
	int y;
	int z;
};

template <class T>
struct component {
	voxel seedVoxel;
	int x;
	int y;
	int z;
	T label;
	T newLabel;
	bool toRelabel;
	int size;
	std::map < T, unsigned int> neighbours;
};

template <class T>
void vtkImageDataCleanLabelsExecute (vtkImageData *input,
									 vtkImageData *output,
									 T *inPtr,
									 T *outPtr)
{
	int dimensions[3];
	input->GetDimensions(dimensions);

	std::vector< component<T> > components;

	int numVoxels=dimensions[0]*dimensions[1]*dimensions[2];

	// set all voxels to not visited (0)
	for (int i=0;i!=numVoxels;i++) {
		outPtr[i]=0;
	}

	int lineOffset=dimensions[0];
	int sliceOffset=dimensions[0]*dimensions[1];

	int maxX=dimensions[0]-1;
	int maxY=dimensions[1]-1;
	int maxZ=dimensions[2]-1;

	std::queue<voxel> voxelsToVisit;

	int x=0,y=0,z=0;
	int neighbourPointer;
	T neighbourLabel;
	T currentLabel;

	// look for connected components
	for (int i=0;i!=numVoxels;i++) {
		if (outPtr[i]==0) {
			currentLabel=inPtr[i];

			// this voxel has not been visited. add it as a new component

			voxel newVoxel;
			newVoxel.pointer=i;
			newVoxel.x=x;
			newVoxel.y=y;
			newVoxel.z=z;

			component<T> newComponent;
			newComponent.label=currentLabel;
			newComponent.seedVoxel=newVoxel;
			newComponent.toRelabel=false;
			int size=0;



			voxelsToVisit.push(newVoxel);
			while (!voxelsToVisit.empty()) {
				voxel currentVoxel=voxelsToVisit.front();
				voxelsToVisit.pop();
				int pointer=currentVoxel.pointer;

				if ((outPtr[pointer]==0)&&(inPtr[pointer]==currentLabel))
				{
					size++;
					int x1=currentVoxel.x;
					int y1=currentVoxel.y;
					int z1=currentVoxel.z;

					voxel neighbourVoxel;
					outPtr[pointer]=1;
					if (x1!=0) {
						neighbourPointer=pointer-1;
						neighbourLabel=inPtr[neighbourPointer];
						if ((outPtr[neighbourPointer]==0)&&(neighbourLabel==currentLabel))
						{
							neighbourVoxel.pointer=neighbourPointer;
							neighbourVoxel.x=x1-1;
							neighbourVoxel.y=y1;
							neighbourVoxel.z=z1;
							voxelsToVisit.push(neighbourVoxel);
						}
						else
						{
							if (neighbourLabel!=currentLabel) {
								newComponent.neighbours[neighbourLabel]++;
							}

						}
					}
					if (x1!=maxX) {
						neighbourPointer=pointer+1;
						neighbourLabel=inPtr[neighbourPointer];
						if ((outPtr[neighbourPointer]==0)&&(neighbourLabel==currentLabel))
						{
							neighbourVoxel.pointer=neighbourPointer;
							neighbourVoxel.x=x1+1;
							neighbourVoxel.y=y1;
							neighbourVoxel.z=z1;
							voxelsToVisit.push(neighbourVoxel);
						}
						else
						{
							if (neighbourLabel!=currentLabel) {
								newComponent.neighbours[neighbourLabel]++;
							}

						}
					}

					if (y1!=0) {
						neighbourPointer=pointer-lineOffset;
						neighbourLabel=inPtr[neighbourPointer];
						if ((outPtr[neighbourPointer]==0)&&(neighbourLabel==currentLabel)) {
							neighbourVoxel.pointer=neighbourPointer;
							neighbourVoxel.x=x1;
							neighbourVoxel.y=y1-1;
							neighbourVoxel.z=z1;
							voxelsToVisit.push(neighbourVoxel);
						}
						else {
							if (neighbourLabel!=currentLabel) {
								newComponent.neighbours[neighbourLabel]++;
							}

						}
					}
					if (y1!=maxY) {
						neighbourPointer=pointer+lineOffset;
						neighbourLabel=inPtr[neighbourPointer];
						if ((outPtr[neighbourPointer]==0)&&(neighbourLabel==currentLabel)) {
							neighbourVoxel.pointer=neighbourPointer;
							neighbourVoxel.x=x1;
							neighbourVoxel.y=y1+1;
							neighbourVoxel.z=z1;
							voxelsToVisit.push(neighbourVoxel);
						}
						else {
							if (neighbourLabel!=currentLabel) {
								newComponent.neighbours[neighbourLabel]++;
							}

						}
					}

					if (z1!=0) {
						neighbourPointer=pointer-sliceOffset;
						neighbourLabel=inPtr[neighbourPointer];
						if ((outPtr[neighbourPointer]==0)&&(neighbourLabel==currentLabel)) {
							neighbourVoxel.pointer=neighbourPointer;
							neighbourVoxel.x=x1;
							neighbourVoxel.y=y1;
							neighbourVoxel.z=z1-1;
							voxelsToVisit.push(neighbourVoxel);
						}
						else {
							if (neighbourLabel!=currentLabel) {
								newComponent.neighbours[neighbourLabel]++;
							}
						}
					}
					if (z1!=maxZ) {
						neighbourPointer=pointer+sliceOffset;
						neighbourLabel=inPtr[neighbourPointer];
						if ((outPtr[neighbourPointer]==0)&&(neighbourLabel==currentLabel)) {
							neighbourVoxel.pointer=neighbourPointer;
							neighbourVoxel.x=x1;
							neighbourVoxel.y=y1;
							neighbourVoxel.z=z1+1;
							voxelsToVisit.push(neighbourVoxel);
						}
						else {
							if (neighbourLabel!=currentLabel) {
								newComponent.neighbours[neighbourLabel]++;
							}

						}
					}
				}
			}

			newComponent.size=size;
			components.push_back(newComponent);
		}

		x++;
		if (x==lineOffset) {
			x=0;
			y++;
			if (y==dimensions[1]) {
				y=0;
				z++;
			}
		}
	}
	cout<<"number of components : "<<components.size()<<endl;
/*	for (unsigned int i=0;i<components.size();i++) {
		cout<<"Component : "<<i<<" : label="<<(int) components[i].label <<", size="<<components[i].size<<endl;
		typename std::map<T , unsigned int>::iterator it;
		cout<<"   neighbours : ";
		for ( it=components[i].neighbours.begin() ; it != components[i].neighbours.end(); it++ )
		{
			cout<<(int) (*it).first<<"("<<(int) (*it).second<<") ";
		}
		cout<<endl;
	}*/

	// find the max label id
	T maxLabel=components[0].label;
	for ( unsigned int i = 0 ; i != components.size() ; i++) {
		T label = components[i].label;
		if ( label > maxLabel ) {
			maxLabel = label;
		}
	}

	cout<<"max label : "<<(int) maxLabel<<endl;

	// fill a vector of labeled components
	std::vector < std::vector < int > > labels;
	labels.resize( maxLabel +1);
	for (unsigned int i = 0 ; i != components.size() ; i++) {
		unsigned int label = static_cast < unsigned int > ( components[i].label );
		labels[ label ].push_back( i );
	}

	std::queue<unsigned int> componentsToRemove;
	for (unsigned int i=0 ; i != labels.size() ; i++) {
		std::vector< int > currentLabelComponents=labels[i];
		if (currentLabelComponents.size()>1) {
			cout<<"label "<<i<<" has "<<currentLabelComponents.size()<<" components"<<endl;

			// find the biggest component
			unsigned int biggestComponent=0;
			int maxSize=0;
			for (unsigned int j=0;j!=currentLabelComponents.size();j++) {
				unsigned int component=currentLabelComponents[j];
				int size=components[component].size;
				if (size>maxSize) {
					biggestComponent=component;
					maxSize=size;
				}
			}
			cout<<"biggest component : "<<biggestComponent<<endl;
			// tag small components to delete
			for (unsigned int j=0;j!=currentLabelComponents.size();j++) {
				unsigned int component=currentLabelComponents[j];
				if (component!=biggestComponent) {
					cout<<"component to delete: "<<component<<endl;
					componentsToRemove.push(component);
				}
			}
		}
	}
	cout<<componentsToRemove.size()<<" component(s) to delete"<<endl;


	while (componentsToRemove.size()) {
		unsigned int componentToRemove=componentsToRemove.front();
		componentsToRemove.pop();
		std::map < T, unsigned int> neighbours=components[componentToRemove].neighbours;
		typename std::map<T , unsigned int>::iterator it;

		// find neighbour with biggest common boundary;
		T neighbourLabel,bestNeighbourLabel=0;
		int maxBoundarySize=-1, boundarySize;
		for ( it=neighbours.begin() ; it != neighbours.end(); it++ )
		{
			neighbourLabel= (*it).first;
			boundarySize= (*it).second;
			if (boundarySize>maxBoundarySize) {
				bestNeighbourLabel=neighbourLabel;
				maxBoundarySize=boundarySize;
			}
		}
		components[componentToRemove].newLabel=bestNeighbourLabel;
		components[componentToRemove].toRelabel=true;
	}


	// find an unused label to use as "not visited" tag
	T unusedLabel=0;
	bool found=false;
	for (unsigned int i=0;i!=(unsigned int) maxLabel+1;i++) {
		if (labels[ i ].size()==0) {
			unusedLabel=i;
			found=true;
			break;
		}
	}
	if (!found) {
		unusedLabel=maxLabel+1;
	}

	// clear output
	// set all voxels to not visited (unusedLabel)
	for (int i=0;i!=numVoxels;i++) {
		outPtr[i]=unusedLabel;
	}

	T paintLabel;

	// fill output with correct labels
	for (unsigned int i=0;i!=components.size();i++) {
		currentLabel=components[i].label;
		if (components[i].toRelabel) {
			paintLabel=components[i].newLabel;
		}
		else {
			paintLabel=components[i].label;
		}

		voxel newVoxel=components[i].seedVoxel;

		voxelsToVisit.push(newVoxel);
		while (!voxelsToVisit.empty()) {
			voxel currentVoxel=voxelsToVisit.front();
			voxelsToVisit.pop();
			int pointer=currentVoxel.pointer;

			if ((outPtr[pointer]==unusedLabel)&&(inPtr[pointer]==currentLabel))
			{
//				size++;
				int x1=currentVoxel.x;
				int y1=currentVoxel.y;
				int z1=currentVoxel.z;

				voxel neighbourVoxel;
				outPtr[pointer]=paintLabel;
				if (x1!=0) {
					neighbourPointer=pointer-1;
					neighbourLabel=inPtr[neighbourPointer];
					if ((outPtr[neighbourPointer]==unusedLabel)&&(neighbourLabel==currentLabel))
					{
						neighbourVoxel.pointer=neighbourPointer;
						neighbourVoxel.x=x1-1;
						neighbourVoxel.y=y1;
						neighbourVoxel.z=z1;
						voxelsToVisit.push(neighbourVoxel);
					}
				}
				if (x1!=maxX) {
					neighbourPointer=pointer+1;
					neighbourLabel=inPtr[neighbourPointer];
					if ((outPtr[neighbourPointer]==unusedLabel)&&(neighbourLabel==currentLabel))
					{
						neighbourVoxel.pointer=neighbourPointer;
						neighbourVoxel.x=x1+1;
						neighbourVoxel.y=y1;
						neighbourVoxel.z=z1;
						voxelsToVisit.push(neighbourVoxel);
					}
				}

				if (y1!=0) {
					neighbourPointer=pointer-lineOffset;
					neighbourLabel=inPtr[neighbourPointer];
					if ((outPtr[neighbourPointer]==unusedLabel)&&(neighbourLabel==currentLabel)) {
						neighbourVoxel.pointer=neighbourPointer;
						neighbourVoxel.x=x1;
						neighbourVoxel.y=y1-1;
						neighbourVoxel.z=z1;
						voxelsToVisit.push(neighbourVoxel);
					}
				}
				if (y1!=maxY) {
					neighbourPointer=pointer+lineOffset;
					neighbourLabel=inPtr[neighbourPointer];
					if ((outPtr[neighbourPointer]==unusedLabel)&&(neighbourLabel==currentLabel)) {
						neighbourVoxel.pointer=neighbourPointer;
						neighbourVoxel.x=x1;
						neighbourVoxel.y=y1+1;
						neighbourVoxel.z=z1;
						voxelsToVisit.push(neighbourVoxel);
					}
				}

				if (z1!=0) {
					neighbourPointer=pointer-sliceOffset;
					neighbourLabel=inPtr[neighbourPointer];
					if ((outPtr[neighbourPointer]==unusedLabel)&&(neighbourLabel==currentLabel)) {
						neighbourVoxel.pointer=neighbourPointer;
						neighbourVoxel.x=x1;
						neighbourVoxel.y=y1;
						neighbourVoxel.z=z1-1;
						voxelsToVisit.push(neighbourVoxel);
					}
				}
				if (z1!=maxZ) {
					neighbourPointer=pointer+sliceOffset;
					neighbourLabel=inPtr[neighbourPointer];
					if ((outPtr[neighbourPointer]==unusedLabel)&&(neighbourLabel==currentLabel)) {
						neighbourVoxel.pointer=neighbourPointer;
						neighbourVoxel.x=x1;
						neighbourVoxel.y=y1;
						neighbourVoxel.z=z1+1;
						voxelsToVisit.push(neighbourVoxel);
					}
				}
			}
		}
	}
}


void vtkImageDataCleanLabels::SimpleExecute (vtkImageData *input, vtkImageData *output)
{
	void* inPtr = input->GetScalarPointer();
	void* outPtr = output->GetScalarPointer();
	switch (input->GetScalarType())
	{
		vtkTemplateMacro(vtkImageDataCleanLabelsExecute(input, output,
										(VTK_TT*)(inPtr), (VTK_TT*)(outPtr)));
	default:
		vtkGenericWarningMacro("ThreadedRequestData: Unknown output ScalarType");
		return;
	}
}

