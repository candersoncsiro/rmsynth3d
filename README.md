# rmsynth3d
Repository for code that generates and visualizes pseudo-3D FITS images from radio astronomy polarization data. Its initial release and comments on its usage are in Rudnick et al. (2024, MNRAS, submitted).

The repository currently consists of two scripts:

1.	generate_pseudo3D.py
   
INPUTS: two matched FITS images in RA, Dec, of polarized intensity and Faraday depth.   These images can be produced either by 
a.	mapping the peak intensity in a Faraday spectrum, and the corresponding Faraday depth at that peak, or
b.	fitting the polarization angle as a function of wavelength-squared to derive a rotation measure and determining the polarized intensity at some fiducial wavelength at each RA,Dec.

OUTPUT: A cube in RA,Dec,Faraday depth space, with values corresponding to the polarized intensity in each voxel.
                
2.	 visualize_pseudo3D.py
   
INPUT:  A cube in RA,Dec,Faraday depth space, with values corresponding to the polarized intensity in each voxel.   This can be produced either through the script generate_pseudo3D.py or from Faraday synthesis.

OUTPUT:  An interactive volumetric rendering of the cube which can be viewed from different angles and with various visualization schemes.

The codes have undergone initial vetting for a variety of inputs, and with a variety of python installations.   However, it will depend on community use and feedback to a) improve and insure robustness for inputs from different telescopes and software packages, and b) adjust default and variable input parameters and create additional analysis and visualization tools.   The available visualization functionality depends on details of the installation for the mayavi package.  Feedback on any issues or ideas for further development should be sent to craig.anderson@anu.edu.au>
