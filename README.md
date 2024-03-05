# rmsynth3d
Repository for code that generates and visualizes pseudo-3D FITS images from radio astronomy polarization data. Its initial release and comments on its usage are in Rudnick et al. (2024, MNRAS, submitted).

The repository currently consists of two scripts, written in python 3:

1.	generate_pseudo3D.py
   
INPUTS: two matched FITS images in RA, Dec, of polarized intensity and Faraday depth.   These images can be produced either by 
a.	mapping the peak intensity in a Faraday spectrum, and the corresponding Faraday depth at that peak, or
b.	fitting the polarization angle as a function of wavelength-squared to derive a rotation measure and determining the polarized intensity at some fiducial wavelength at each RA,Dec.

OUTPUT: A cube in RA,Dec,Faraday depth space, with values corresponding to the polarized intensity in each voxel.
                
2.	 visualize_pseudo3D.py
   
INPUT:  A cube in RA,Dec,Faraday depth space, with values corresponding to the polarized intensity in each voxel.   This can be produced either through the script generate_pseudo3D.py or from Faraday synthesis.

OUTPUT:  An interactive volumetric rendering of the cube which can be viewed from different angles and with various visualization schemes.

The codes have undergone initial vetting for a variety of inputs, and with a variety of python installations.   However, it will depend on community use and feedback to a) improve and insure robustness for inputs from different telescopes and software packages, and b) adjust default and variable input parameters and create additional analysis and visualization tools.   The available visualization functionality depends on details of the installation for the mayavi package.  Feedback on any issues or ideas for further development should be sent to craig.anderson@anu.edu.au>

___________________________________________________________________________________________________________________________

Docstring for generate_pseudo3D.py:

"""
generate_pseudo3D.py. Written by Craig Anderson and Larry Rudnick in Feb 2024 and described in Rudnick et al., 2024, MNRAS, submitted.

Generates a FITS cube from FITS maps of polarized intensity and Faraday depth, with user-specified Gaussian smoothing in Faraday depth.  The output FITS cube can be viewed in any appropriate software package, such as ds9 or CARTA, or used as input to the visualize_pseudo3D.py script.

Usage:
	python generate_pseudo3D.py p_filename.fits rm_filename.fits output_cube.fits 
--rm_min -100 --rm_max 100 --rm_incr 1 --smooth_sigma 3

Parameters:
	p_filename (str): Path to the polarized intensity FITS file.
	rm_filename (str): Path to the Faraday depth FITS file.
	output_filename (str): Path for the output FITS cube file.
	rm_min (int, rad/m^2): Minimum Faraday depth value.
	rm_max (int, rad/m^2): Maximum Faraday depth value.
	rm_incr (float, rad/m^2): Pixel size in Faraday depth space.
	smooth_sigma (odd integer, channels): Smoothing factor (standard deviation of the Gaussian kernel) applied along  the Faraday depth axis for the pseudo 3D representation.

The file names have no defaults.
					  
NOTE: This script has been tested and verified on FITS maps that either have NAXIS=2, or NAXIS=3 with a singleton degenerate 'LINEAR' type axis (shape: <4000,<4000,1). Development work to handle other possible cases is most welcome!

Installation requirements:  
numpy, astropy, scipy
"""

___________________________________________________________________________________________________________________________

Docstring for visualize_pseudo3D.py:

"""
Generate Pseudo3D Visualizations from FITS Cubes.

This script, generate_pseudo3D.py, produces 3D volumetric renders of FITS cubes of RA, Decl, Faraday depth. Developed by Craig Anderson and Larry Rudnick, February 2024. This script is a companion to 'generate_pseudo3D.py' and supports the findings presented in Rudnick+2024.

Usage:
    python visualize_pseudo3D.py <file_path> --rsf X Y Z

Example:
    python visualize_pseudo3D.py /path/to/data_cube.fits --rsf 0.5 0.5 2

Parameters:
    file_path (str): Path to the FITS cube file.
    --rsf (X Y Z): Resampling factors for RA, Dec, and Faraday depth, respectively.

Note: For RA and Dec, only resampling factors <1 are allowed, and will reduce the number of pixels along those axes by that multiplicative factor, with the size of each pixel increased.
For Faraday depth, factors may be <1 or >1. For an rsf>1, the number of pixels along the Faraday depth axis will be increased by that multiplicative value, and the size of each pixel correspondingly decreased.
Choices of resampling factors may be useful to manage computational resources and to improve the visualization of structures in Faraday space. 

Features:
- Utilizes Mayavi for interactive 3D visualization, offering up to the full feature set of this package, depending on installation. Mayavi has known installation issues that can render some advanced functionality inoperable. When available, a robust recipe for obtaining the full mayavi functionality will be provided here.

Installation requirements:  
	numpy, astropy, scipy, mayavi (which requires vtk and a GUI toolkit)*

    
*See, e.g.,  https://docs.enthought.com/mayavi/mayavi/installation.html or https://github.com/prabhuramachandran/mayavi-tutorial/blob/master/installation.md
"""
