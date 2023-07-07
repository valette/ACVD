import nibabel as nib
import os

def flipAndSaveToRAS(filename, args):
    """
    Function used to flip the orientation of the original data to RAS+ orientation
    Save the data flipped under the same name into a sub-directory called "RASData"
    @input filename : path of the file to flip to RAS
    """
    
    #Recover the image object
    img = nib.load(filename)
    
    #Get the current orientation
    CurrentOrientation = nib.aff2axcodes(img.affine)
    print("The current orientation is : ", CurrentOrientation)
    
    #Check if the current orientation is already RAS+
    if CurrentOrientation == ('R', 'A', 'S') :
        
        print("Image already recorded into the RAS+ orientation, nothing to do")
        
    else :
        #Flip the image to RAS
        img = nib.as_closest_canonical(img)
                
    ##Check the new orientation
    NewOrientation = nib.aff2axcodes(img.affine)
        
    ##Save the image into the folder RASData with the same name
    index = filename.rfind('/')
    outputDirectory = os.path.join(filename[0:index+1],'RASData');
    if args.outputDirectory : outputDirectory = args.outputDirectory 

    if not os.path.exists( outputDirectory ): os.mkdir( outputDirectory )
        
    #Set Qcode to 1 that the Qform matrix can be used into the further processing
    img.header['qform_code'] = 1
        
    #Save the flipped image
    outputFileName = filename[index+1:]
    if args.outputFileName : outputFileName = args.outputFileName
    fullFileName = os.path.join(outputDirectory, outputFileName )
    img_data = img.get_fdata()
    img_conv = nib.Nifti1Image(img_data.astype(img.header.get_data_dtype()), img.affine, img.header)
    nib.save(img_conv, fullFileName)
        
    print("The new orientation is now : ", NewOrientation)
        
    #########
    ## Test #
    #########
        
    ###Check if the we saved the RAS+ data
    #imTest = nib.load(fullFileName)
    #print(nib.aff2axcodes(imTest.affine), imTest.get_qform(), imTest.header['qform_code'])
                

def getProgramParameters():
    """
    Function used to get the filename passed as parameters while calling the program
    @output args.filename : path passed as parameter
    """
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Flip the orientation of the file (-f) or directory (-d) to RAS+ and save it into a sub-directory called RASData')
    parser.add_argument('-f', dest = 'filename', help = 'File which orientation is to change')
    parser.add_argument('-d', dest = 'directory', help = 'Directory where the files are located')
    parser.add_argument('-o', dest = 'outputDirectory', help = 'Output directory where the files will be written')
    parser.add_argument('-of', dest = 'outputFileName', help = 'Output file name')
    args = parser.parse_args()
    
    return args
    
def main():
    """
    Function to run as main routine
    """
    #Recover the parameters passed as entry of the program
    args = getProgramParameters()
    #Recover the files to flip
    if args.directory :
        listFiles = [ os.path.join(args.directory, f) for f in os.listdir(args.directory) if f.endswith(".nii.gz")]
    elif args.filename :
        listFiles = [args.filename]
    else : 
        print("Need to specify a file (-f) or a directory (-d) \n Processus aborted")
        return 0
    
    for filename in listFiles :
        print("\nFlip ", filename)
        #Call to the function that flip the input data to RAS model and save it
        flipAndSaveToRAS(filename, args)
    
    
if __name__ == '__main__':
    main()
