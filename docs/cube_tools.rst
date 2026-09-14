Cube Utility Tools
==================

The :mod:`CubeGen.cubetools.cube_tools` module provides utility
functions for analysing and validating reconstructed data cubes.

Astrometric Matching
--------------------

``astromatch`` compares the astrometric registration of two reconstructed
data cubes using the centroid of a source measured from an integrated
flux image.

For each cube, the spectral axis is collapsed to construct a
two-dimensional image:

.. math::

   I(x,y) = \sum_{\lambda} F(\lambda,x,y).

A two-dimensional PSF is fitted to the resulting image in order to
determine the source centroid.

The centroid measured in the reference cube is then converted to
celestial coordinates using its WCS and projected onto the pixel
coordinate system of the second cube.

The residual astrometric displacement is

.. math::

   \Delta x_{\rm WCS} = x_1 - x_{0\rightarrow1}

and

.. math::

   \Delta y_{\rm WCS} = y_1 - y_{0\rightarrow1},

where :math:`x_{0\rightarrow1}` and :math:`y_{0\rightarrow1}` are the
expected coordinates of the reference source in the second cube.

Values close to zero indicate that the measured source positions are
consistent with the WCS registration of the cubes.


Usage
~~~~~

.. code-block:: python

   from CubeGen.cubetools.cube_tools import astromatch

   dx, dy, dx_wcs, dy_wcs = astromatch(
       "cube_reference.fits",
       "cube_target.fits",
       sig=2
   )

   print("Direct pixel offset:", dx, dy)
   print("WCS residual:", dx_wcs, dy_wcs)


API Reference
-------------

.. autofunction:: CubeGen.cubetools.cube_tools.astromatch