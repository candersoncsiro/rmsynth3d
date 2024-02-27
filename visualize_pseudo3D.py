"""
`visualize_pseudo3D.py`: A tool for 3D volumetric visualization of RA, Decl, RM FITS data cubes. Developed by Craig Anderson and Larry Rudnick, February 2024. This script is a companion to 'generate_pseudo3D.py' and supports the findings presented in Rudnick+2024.

Usage example: `python visualize_pseudo3D.py path/to/your/data_cube.fits --rsf 1.5 0.25 0.25`

Positional arguments:
- `file_path` (str): The path to the FITS file intended for 3D rendering.

Optional arguments:
--rsf X Y Z allows for specifying multiplicative resampling factors X Y Z for each dimension (RM, RA, Decl) of the data cube, adjusting the number of pixels in each dimension to refine the visualization or manage computational resources. For instance, `--rsf 1.5 0.25 0.25` will increase the RM dimension by a factor of 1.5, while reducing the RA and Dec dimensions by a factor of 4, enhancing certain visualization aspects or alleviating memory constraints. Note: Resampling factors above 1 for RA and Dec are discouraged due to potential for long runtimes and memory issues.

Features:
- Utilizes Mayavi for interactive 3D visualization, offering the full feature set of this package.
"""

import numpy as np
import argparse
from astropy.io import fits
from mayavi import mlab
from scipy.ndimage import zoom

def resample_cube(cube, factors):
	"""
	resample the 3D data cube by the specified factor for each dimension.
	Parameters:
	cube (ndarray): The original 3D numpy array to be resampled.
	factors (tuple): A tuple of three numbers indicating the multiplicative resampling factor for each dimension.
	Returns:
	ndarray: The resampled 3D numpy array.
	"""
	# The following 'if' statement, and forced resampling (even with factors (1, 1, 1)  is a kludge to to avoid a strange bug in in VTK if the cube is not re-sampled, even for unit re-sample factors (i.e. no nominal change to cube). See also 'load_fits_data' function. Moving this out of "if" statement as a kludge fix. (Error is: site-packages/tvtk/util/ctf.py:123: RuntimeWarning: invalid value encountered in double_scalars). 
	if any(f != 1 for f in factors):
		print(f"\nResampling cube by factors {factors} along RA, Dec, RM axes...")
	else:
		pass
	return zoom(cube, zoom=factors, order=0)  # Using order=0 for NN interpolation

def validate_resample_factors(factors):
	"""
	Validates that the resampling factors meet the specified requirements.
	"""
	if any(f > 1 for f in factors[:-1]):
		raise ValueError("Invalid resampling factors. RA and Dec factors must be <= 1.")

def check_for_nans_infs(cube):
	"""
	Check the data cube for NaNs and infinite values.

	Parameters:
	cube (ndarray): A 3D numpy array representing the data cube.

	Raises:
	ValueError: If NaNs or infinite values are found in the data cube.
	"""
	print('\nChecking for NaNs, infs...')
	if np.isnan(cube).any() or np.isinf(cube).any():
		raise ValueError("Data cube contains NaNs or infinite values. Aborting.")    	

def render_3d_cube(cube):
	"""
	Render a semi-transparent 3D volumetric representation of the 3D data cube.
	Parameters:
	cube (ndarray): A 3D numpy array representing the data cube to be rendered.
	"""
	check_for_nans_infs(cube)  # Check for NaNs and Infs before rendering
	print(f"Data range of cube: min={cube.min():.2g}, max={cube.max():.2g}")
	print('Rendering cube output...')
	src = mlab.pipeline.scalar_field(cube)
	vol = mlab.pipeline.volume(src, vmin=cube.min(), vmax=cube.max())
	vol._volume_property.scalar_opacity_unit_distance = 0.1
	#View orientation init --- looking down the FD axis (neg->pos)
	mlab.view(azimuth=0, elevation=0)  
	mlab.roll(0)
	mlab.orientation_axes()
	mlab.show()

def calculate_resample_factors(original_dims, max_pixels=2000):
	"""
	Calculate the resample factors needed to ensure dimensions do not exceed max_pixels.

	Parameters:
	original_dims (tuple): The original dimensions of the data cube (RM, RA, Dec).
	max_pixels (int): The maximum allowed pixels in the RA and Dec dimensions.

	Returns:
	tuple: Calculated resample factors for each dimension to meet the max_pixels constraint.
	"""
	factors = []
	for dim in original_dims:
		if dim > max_pixels:
			factor = max_pixels / dim
		else:
			factor = 1
		factors.append(factor)
	return tuple(factors)

def load_fits_data(file_path, user_resample_factors):
	"""
	Load the 3D data cube from a FITS file and check if the specified resample factors will keep RA and Dec dimensions within the desired limits.
	Throws an exception if the resampling factors do not meet the constraints.

	Parameters:
	file_path (str): Path to the FITS file containing the 3D data cube.
	user_resample_factors (tuple): User-specified resampling factors for each dimension.

	Returns:
	ndarray: A 3D numpy array loaded from the specified FITS file.
	"""
	print('Loading FITS...')
	with fits.open(file_path) as hdul:
		data = hdul[0].data
		# deal with Astropy / Astropy FTIS axis ordering convention mismatch
		data_transposed = np.transpose(data, (2, 1, 0)) 
		original_dims = data_transposed.shape  # Original dimensions of the data cube (RM, RA, Dec)

	# Calculate necessary resample factors to keep dimensions within limits that are less likely to produce OOM errors, long runtimes, etc.
	necessary_factors = calculate_resample_factors(original_dims[1:], max_pixels=2000)  # Skip RM dimension for this calculation

	# Check if user's resample factors are sufficient
	for i, factor in enumerate(user_resample_factors[1:], start=1):  # Again, skip RM dimension
		if original_dims[i] * factor > 2000:
			raise ValueError(f"Specified or default resample factor for dimension {i} ({factor}) will result in more than 2000 pixels. "
							 "Please adjust the --rsf flag to ensure dimensions do not exceed this limit. "
							 f"Recommended maximum factor for this dimension: {necessary_factors[i-1]:.2f}.")

	# If this point is reached, the user's factors are deemed sufficient, or no resampling is needed
	print(user_resample_factors)
	if any(f != 1 for f in user_resample_factors):
		print("User-specified resampling factors are within acceptable limits.")
	
	#Currently, there is a strange bug that triggers deep in VTK if the cube is not re-sampled, even for unit re-sample factors (i.e. no nominal change to cube). Moving this out of "if" statement as a kludge fix. (Error is: site-packages/tvtk/util/ctf.py:123: RuntimeWarning: invalid value encountered in double_scalars)
	data_transposed = resample_cube(data_transposed, user_resample_factors)
	
	return data_transposed

def main():
	"""
	Main function to parse command line arguments and render the 3D data cube from a FITS file.
	"""
	parser = argparse.ArgumentParser(description=globals()['__doc__'])
	parser.add_argument("file_path", help="Path to the FITS file containing the 3D data cube")
	parser.add_argument("--rsf", nargs=3, type=float, default=[1, 1, 1],
						help="Multiplicative resampling factors for each dimension (RA, Decl), e.g., 0.2 0.2 2 will reduce the number of pixels in RA and Dec by a factor of 5 in each case, and increase the number of pixels in the RM dim by a factor of 2 (expanding it in the viz cube). Factors above 1 are not alowed for the RA and Dec dims, becuase they tend to produce long runtimes and Out Of Memory (OOM) errors.")
	args = parser.parse_args()

	# Validate rsf before proceeding
	validate_resample_factors(args.rsf)  

	# Load and optionally resample the data before processing further
	cube = load_fits_data(args.file_path, tuple(args.rsf))
	render_3d_cube(cube)

if __name__ == "__main__":
	main()
