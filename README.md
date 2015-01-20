# roi-extraction
DESCRIPTION

This toolbox allows a semi-automated roi definition procedure.
The graphics user interface requires that you input:

1 - the folder in which your subjects are located (e.g. /home/myproject/mridata)

2 - the SPM.mat of the first subject for the functional localiser from which you want to expand ROIs (e.g. /home/myproject/mridata/subject1/functionallocaliserstats/SPM.mat)

3 - the native space anatomical of the first subject of the study (e.g. /home/myproject/mridata/subject1/structurals/sMPRAGE.nii)

4 - an anatomical mask in which the search will be restrained.

5 - whether you want to use contiguous expansion of the ROIs or discontiguous 
(discontiguous is not recommended for now and an option to add an anatomically defined inclusion mask will be ported in future releases)

6 - the Number of voxels that you wish to save in your binary mask (this also handles matlab regular expression e.g. 
linspace(20,400,20) or vectors e.g. [20 40 60 80 100])

7 - the activation contrast that you entered in SPM at the first level to define this ROI

8 - A name that you want to use for your ROI masks (e.g. FFA). The name will be used in twofolds:

a: a folder will be created to store binary .img masks with this name (e.g. /home/myproject/mridata/subject1/masks/FFA/)

b: the name will be used as a suffix for the number of voxels chosen (e.g. /home/myproject/mridata/subject1/masks/FFA/leftFFA_20.img,/home/myproject/mridata/subject1/masks/FFA/leftFFA_40.img...)
   
As hinted above, the toolbox will proceed with left and right hemisphere definitions independently, one after the other. The toolbox will recursively select the subjects contained in your mridata folder so it is important that your folder architecture remain the same across all the subjects.
    
REQUIREMENTS:

   - SPM 8 OR ABOVE NEEDS TO BE DEFINED IN MATLAB PATH

Ian Charest - August 2014

Medical Research Council - Cognition and Brain Sciences Unit

 
Last Modified by Ian Charest 20-Jan-2015 11h17 am
