Reconstruction Algorithm
========================

3DCubeGen reconstructs regular 3D data cubes (RA, Dec, wavelength) from
irregularly sampled IFU fiber data using weighted kernel interpolation. The
same core algorithm is shared across all three supported spectrographs (LVM,
MEGARA, VirusP).


Kernel Interpolation
--------------------

For each output spaxel, the algorithm computes a weighted sum of all fiber
spectra within a cutoff radius. The weight assigned to each fiber depends on
its distance from the spaxel center:

.. math::

   w_k = \exp\!\left(-\frac{1}{2}\left(\frac{R_k}{\sigma}\right)^{\!\alpha}\right)

where:

- :math:`R_k` is the Euclidean distance from fiber *k* to the spaxel center
  (arcsec)
- :math:`\sigma` is the smoothing kernel size (``sigm_s`` config parameter,
  arcsec)
- :math:`\alpha` is the shape factor (``alph_s`` config parameter)

The shape factor controls the kernel profile:

- :math:`\alpha = 2` produces a **Gaussian** kernel (the default)
- :math:`\alpha > 2` produces a sharper, more top-hat-like kernel

The interpolated flux at each spaxel is the weighted mean:

.. math::

   F_{\text{spaxel}} = \frac{\sum_k F_k \, w_k}{\sum_k w_k}

Only fibers with finite, positive flux values contribute to the sum. Spaxels
with zero total weight are left at zero flux.


Radius Cutoff
-------------

To avoid unnecessary computation, only fibers within a cutoff radius of each
spaxel are considered. The cutoff radius is:

.. math::

   r_{\text{cut}} = \max\!\left(\frac{\sigma}{2},\; \frac{d_{\text{fib}} \times 3.5 \times 2}{2}\right)

where :math:`d_{\text{fib}}` is the fiber diameter. For LVM, the fiber
diameter parameter (``fibA``) is 35.3 arcsec; for VirusP it is configurable
(default 4.16 arcsec); for MEGARA it is derived from the hexagonal bundle
geometry.

Additionally, a Y-axis pre-filter selects only fibers whose Y coordinate is
within :math:`r_{\text{cut}}` of the current row, reducing the inner loop
to nearby fibers.


Error Propagation
-----------------

When error propagation is enabled (``errors: True``), errors are propagated
using inverse-variance weighting:

.. math::

   \sigma_{\text{spaxel}} = \frac{\sqrt{\sum_k (\sigma_k \, w_k)^2}}{\sum_k w_k}

where :math:`\sigma_k` is the error spectrum of fiber *k*. This preserves the
correct noise scaling through the weighted interpolation.


ADR-Aware Mode
--------------

When fiber positions vary as a function of wavelength (e.g., due to
Atmospheric Differential Refraction correction), the position arrays become
2D: ``x_ifu_V[n_fibers, n_wavelengths]``. The algorithm detects this
automatically (``len(x_ifu_V.shape) > 1``) and switches to per-wavelength
distance calculations:

- In **1D mode** (fixed positions): a single scalar distance is computed per
  fiber-spaxel pair, and the weight applies uniformly across the spectrum.
- In **2D mode** (wavelength-dependent positions): the distance
  :math:`R_k(\lambda)` is computed for each wavelength channel independently,
  producing wavelength-dependent weights.

This mode is used by the VirusP pipeline when ``adrcor: True``.


Spaxel Size
-----------

The output spaxel size is controlled by the ``pix_s`` parameter. When set to
``0`` (auto-calculate), each spectrograph uses its own default:

- **LVM**: ``pix_s = 0.75 * sigm_s``
- **MEGARA**: ``pix_s = sigm_s``
- **VirusP**: ``pix_s = 0.75 * sigm_s``


Flux Scaling
------------

After interpolation, the flux in each spaxel is scaled by an area factor to
convert from fiber-integrated flux to surface brightness:

.. math::

   \text{factor} = \frac{(\text{pix\_s})^2}{\pi \, (d_{\text{fib}} / 2)^2}

This accounts for the different solid angles subtended by a square spaxel
versus a circular fiber aperture.


Multi-Threading
---------------

The reconstruction is parallelized along the spatial Y-axis using Python's
``ThreadPool``. The number of threads is controlled by the ``nproc`` config
parameter. Each thread processes a contiguous block of rows independently.

The outer loop iterates over X columns, and for each column the Y-axis work is
split across ``nproc`` threads, with each thread calling the ``ifu_const()``
kernel function on its assigned row range.


Implementation Reference
------------------------

The core interpolation functions are in ``CubeGen/tools/kernel.py``:

- ``ifu_const()`` -- cube reconstruction kernel (flux + optional error)
- ``map_const()`` -- map reconstruction kernel (flux + magnitude + errors)
- ``matrix_const()`` -- matrix-mode variant used for covariance estimation

The kernel functions are called by the spectrograph-specific drivers:

- ``CubeGen/tools/cubegen.py::map_ifu()`` -- LVM cube generation
- ``CubeGen/tools/mapgen.py::gen_map()`` -- LVM map generation
- ``CubeGen/megaratools/megcubegen.py::megmap_ifu()`` -- MEGARA cube generation
- ``CubeGen/virustools/viruscubegen.py::virusmap_ifu()`` -- VirusP cube generation
