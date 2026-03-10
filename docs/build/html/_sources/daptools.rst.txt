DAP Integration Tools
=====================

The ``CubeGen.daptools`` module provides utilities for interfacing
reconstructed cubes with downstream analysis tools, particularly the
Data Analysis Pipeline (DAP) used by the LVM survey.


RSS Extraction from Cubes
--------------------------

**Module:** ``CubeGen.daptools.rss_gen``

The ``rssp_extract()`` function converts a reconstructed 3D cube back into a
pseudo-RSS (Row-Stacked Spectra) format. This is used when the downstream
analysis pipeline expects RSS-format input rather than 3D cubes.

.. code-block:: python

   from CubeGen.daptools.rss_gen import rssp_extract

   rssp_extract(
       name='NGC_1234',
       path='./cubes/',
       path_out='./rss/',
       basename_in='lvmCube-NAME.fits.gz',
       basename_out='lvmRSS-NAMElab.fits',
       flu16=True,
       nsplit=0,
   )

**How it works:**

1. Reads the input cube (flux + error extensions)
2. Iterates over all spatial pixels (spaxels)
3. Each non-zero spaxel becomes one "fiber" in the output RSS
4. Constructs a SLITMAP table with RA/Dec coordinates derived from the
   cube's WCS, plus pixel positions for back-referencing
5. Writes the result as a multi-extension FITS file (see :doc:`outputs`)


Producing RSS output from sincube
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``sincube`` command does not have a ``cube2rss`` option -- that feature is
only available in ``pipecube``.  However, ``rssp_extract()`` requires only the
cube file as input, so it can be called manually after ``sincube`` completes:

.. code-block:: python

   from CubeGen.daptools.rss_gen import rssp_extract

   rssp_extract(
       name='NGC_1234',          # matches the name used in the sincube config
       path='out_cubes/',        # out_path from the sincube config
       path_out='out_rss/',
   )

The function reads HDU 0 (flux) and HDU 1 (error) from
``lvmCube-NGC_1234.fits.gz`` and writes ``lvmRSS-NGC_1234.fits.gz`` to
``out_rss/``.

.. important::

   ``sincube`` must be run with ``errors: True`` in the config file.  When
   errors are disabled, HDU 1 of the cube is an array of ones, and the IVAR
   extension of the RSS will be filled with meaningless values
   (``1 / 1.0² = 1`` everywhere).


Which RSS extensions contain real data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``rssp_extract()`` works entirely from the cube file.  Because the cube does
not preserve the original per-fiber sky spectra or LSF, several output
extensions cannot be reconstructed and are filled with placeholder values:

.. list-table::
   :header-rows: 1
   :widths: 10 20 15 55

   * - HDU
     - EXTNAME
     - Status
     - Content
   * - 0
     - ORIGINAL
     - **Real**
     - Primary header copied from the cube, with ``POSCIRA``/``POSCIDE``
       set from the cube's ``CRVAL1``/``CRVAL2``.  ``EXPOSURE`` is
       hardcoded to ``000000``.
   * - 1
     - FLUX
     - **Real**
     - Flux spectra extracted spaxel-by-spaxel from the cube HDU 0.
   * - 2
     - IVAR
     - **Real** (requires ``errors: True``)
     - Inverse variance computed as ``1 / error²`` from cube HDU 1.
       If the cube was generated without error propagation, this
       extension contains ones and should not be used.
   * - 3
     - MASK
     - **Placeholder**
     - Zero array -- no mask information survives cube reconstruction.
   * - 4
     - WAVE
     - **Real**
     - Wavelength array reconstructed from the cube header keywords
       ``CRVAL3``, ``CDELT3``/``CD3_3``, and ``CRPIX3``.
   * - 5
     - LSF
     - **Placeholder**
     - Ones array -- the line spread function is not stored in the cube.
   * - 6
     - SKY_EAST
     - **Placeholder**
     - Ones array -- sky spectra are not preserved in the cube.
   * - 7
     - SKY_EAST_IVAR
     - **Placeholder**
     - Ones array.
   * - 8
     - SKY_WEST
     - **Placeholder**
     - Ones array.
   * - 9
     - SKY_WEST_IVAR
     - **Placeholder**
     - Ones array.
   * - 10
     - SLITMAP
     - **Real**
     - Binary table.  ``ra`` and ``dec`` are computed from the cube WCS
       via ``pixel_to_skycoord()``.  ``xpmm``/``ypmm`` store the pixel
       column/row in the cube for back-referencing.  Bookkeeping columns
       (``targettype``, ``telescope``, ``ifulabel``, etc.) are filled
       with generic placeholder strings (``'science'``, ``'Sci'``,
       ``'Sci1'``).

This is the same behaviour whether the cube was produced by ``sincube`` or
``pipecube`` -- ``rssp_extract()`` is called identically in both cases and
uses no information beyond the cube file itself.


**Spatial tiling support:**

For large cubes that exceed memory limits, the ``nsplit`` parameter enables
spatial tiling. When ``nsplit > 1``, the cube is divided into an
``nsplit x nsplit`` grid, and only the tile specified by ``spt=[i, j]`` is
processed. The output filename is tagged with the tile index (e.g.,
``lvmRSS-NGC_1234_p01.fits``).

**1D RSS from external spectra:**

The ``rssp_multi1d()`` function creates a pseudo-RSS from a collection of
independent 1D spectra (not from a cube). It reads a file list, resamples all
spectra to a common wavelength grid, and writes the result in the same
multi-extension format as ``rssp_extract()``. This is useful for collecting
spectra from different sources into a single DAP-compatible file.


DAP Emission Line Tools
------------------------

**Module:** ``CubeGen.daptools.dap_fluxelines``

These tools extract emission line measurements from DAP analysis products and
map them back onto the 2D spatial grid of the original cube.

``extract_values()``
^^^^^^^^^^^^^^^^^^^^

Reads emission line fit results from the ``NP_ELINES_B``, ``NP_ELINES_R``,
and ``NP_ELINES_I`` extensions of a DAP output file. Extracts for each line:

- Flux and error (``flux``, ``e_flux``)
- Velocity and error (``vel``, ``e_vel``)
- Velocity dispersion and error (``disp``, ``e_disp``)
- Equivalent width and error (``EW``, ``e_EW``)

``extract_PMvalues()``
^^^^^^^^^^^^^^^^^^^^^^

Reads emission line fit results from the ``PM_ELINES`` and ``PM_KEL``
extensions (parametric and kinematic emission line fits). Extracts flux,
velocity, and dispersion for each modeled wavelength.

``get_values()``
^^^^^^^^^^^^^^^^

Combines ``extract_values()`` and ``extract_PMvalues()`` to build complete
2D emission line maps. For each spaxel in the original cube grid, it matches
the DAP fiber ID to the RSS SLITMAP coordinates and places the measurements
at the correct spatial position.

Returns three map arrays:

- ``map_dap`` -- non-parametric line measurements (NP), shape
  ``[n_lines * 8, nx, ny]`` (8 quantities per line)
- ``map_dap2`` -- parametric line measurements (PM), shape
  ``[n_models * 6, nx, ny]``
- ``map_dap3`` -- kinematic emission line measurements (PM_KEL), shape
  ``[n_models * 6, nx, ny]``

``dap_extract_fluxelines()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

High-level wrapper that orchestrates the full extraction workflow:

1. Reads the cube header to determine spatial dimensions and WCS
2. Optionally handles spatial tiling (``nsplit > 1``)
3. Calls ``get_values()`` for each tile or the full cube
4. Writes the output as ``lvmDAPMap-<name>-flux_elines.fits.gz`` with three
   extensions: NP, PM, and PM_KEL

.. code-block:: python

   from CubeGen.daptools.dap_fluxelines import dap_extract_fluxelines

   dap_extract_fluxelines(
       name='NGC_1234',
       out_path='./dap_maps/',
       path_cube='./cubes/',
       path_dap='./dap/',
       path_rss='./rss/',
   )

The output header contains ``VAL_xxxx`` keywords that map each slice of the
data cube to the corresponding emission line and measurement type (e.g.,
``flux_Ha_6565``, ``vel_OIII_5008``).
