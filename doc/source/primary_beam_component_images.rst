Primary beam (component) images
===============================

As described in the :doc:`primary beam correction chapter <primary_beam_correction>`, in most primary-beam correcting modes, WSClean will store 16 separate images for each Mueller-matrix component (e.g. named ``wsclean-beam-0.fits``, ``wsclean-beam-1.fits``, ..., ``wsclean-beam-15.fits``). These are difficult to interpret, as e.g. the shape that any one of these images indicates might not reflect the shape of the beam. 

To turn a Mueller matrix into something that can be interpreted / validated, one option is to use the 16 images to form the Mueller matrix for each pixel (see below) and multiply this with the vector ``(0.5; 0.0; 0.0; 0.5)``. This results in the full Jones sensitivity of the beam. Starting :doc:`WSClean 3.2 <changelogs/v3.2>`, only the components of the Mueller matrix that are necessary for the correction are stored, e.g. with Stokes I imaging, only element 0 and 15 are stored. 

Mueller matrix structure
------------------------
The images form the lower triangle of the Hermitian matrix. For the diagonal elements, only the real value is stored, and off-diagonal elements contain two images: the real and imaginary parts. They're counted from left to right, top to bottom, so their indices correspond to the following positions::

  [0]
  [1]+[ 2]i   [ 3]
  [4]+[ 5]i   [ 6]+[ 7]i   [ 8]
  [9]+[10]i   [11]+[12]i   [13]+[14]i   [15]

This implies for example that the diagonal elements are stored in images named like ``wsclean-beam-0.fits``, ``wsclean-beam-3.fits``, ``wsclean-beam-8.fits``, ``wsclean-beam-15.fits``, and the bottom-left entry of the matrix is split into a real part (``wsclean-beam-9.fits``) and an imaginary part (``wsclean-beam-10.fits``).

