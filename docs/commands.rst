Commands Overview
=================

The ``3dcubegen`` CLI provides four commands for LVM data. They fall into two
categories -- **cube generation** (full 3D spectral cubes) and **map
generation** (2D narrow/broad-band images) -- each available in a single-target
or batch-pipeline mode.

.. code-block:: bash

   3dcubegen <command> -c <config_file>

.. list-table::
   :header-rows: 1
   :widths: 14 14 72

   * - Command
     - Output
     - Purpose
   * - ``sincube``
     - 3D cube
     - Generate one data cube from a single target's exposure list.
   * - ``pipecube``
     - 3D cube(s)
     - Batch-process multiple targets into cubes, with optional spectral
       chunking, co-addition, and RSS extraction.
   * - ``sinmap``
     - 2D map(s)
     - Generate narrow-band or broad-band 2D images for a single target.
   * - ``pipemap``
     - 2D map(s) + color JPEGs
     - Batch-process multiple targets into all seven filter maps and
       composite color images (OHS, OOH, gri, griOOH).


Single vs. Pipeline Commands
-----------------------------

The "sin" (single) and "pipe" (pipeline) variants differ in three key ways:

**1. Number of targets**

- ``sincube`` / ``sinmap`` process **one target** at a time. The target name
  is set by ``nameF``, and its coordinates by ``coord_cen``.
- ``pipecube`` / ``pipemap`` process a **list of targets** in a loop. Target
  names are given as a YAML list in ``nameL``, with matching coordinates in
  ``coord_cenL``.

**2. Kernel specification**

- ``sincube`` / ``sinmap`` accept explicit kernel parameters: ``sigm_s``
  (kernel size in arcsec), ``pix_s`` (spaxel size), and ``alph_s`` (shape
  factor). This gives full control over spatial resolution.
- ``pipecube`` / ``pipemap`` use a single-letter ``type`` code (``a`` through
  ``j``) that maps to predefined ``sigm_s`` and ``pix_s`` values via
  ``kernel_pipe()``. The shape factor ``alph_s`` is fixed at 2.0 (Gaussian).
  See :ref:`kernel-types` for the full table.

**3. Output path structure**

- ``sincube`` / ``sinmap`` write output directly into ``out_path``.
- ``pipecube`` / ``pipemap`` append ``<redux_ver>/<type>/`` to the output
  path automatically (e.g. ``out_cubes/1.1.1.dev0/c/``).


Cube vs. Map Commands
---------------------

**Cube commands** (``sincube``, ``pipecube``) produce full 3D FITS data cubes
(RA x Dec x wavelength). They read RSS frames, interpolate fiber spectra onto
a regular spatial grid using the weighted kernel algorithm, and write the
result as ``lvmCube-<name>.fits.gz``.

**Map commands** (``sinmap``, ``pipemap``) produce 2D FITS images by
convolving the RSS spectra with photometric filter curves before spatial
interpolation. The output files are named ``lvmMap-<name>_<band>.fits.gz``.

The available photometric bands are:

.. list-table::
   :header-rows: 1
   :widths: 8 15 20 57

   * - ID
     - Code
     - Band
     - Type
   * - 0
     - OII
     - [OII] 3728 A
     - Narrow-band (doublet, 8 A width)
   * - 1
     - OIII
     - [OIII] 5008 A
     - Narrow-band (5 A width)
   * - 2
     - HI
     - H-alpha 6565 A
     - Narrow-band (5 A width)
   * - 3
     - SII
     - [SII] 6733 A
     - Narrow-band (5 A width)
   * - 4
     - g
     - SDSS g
     - Broad-band
   * - 5
     - r
     - SDSS r
     - Broad-band
   * - 6
     - i
     - SDSS i
     - Broad-band

Map-specific parameters:

- ``redshift`` -- target redshift, used to shift filter curves to rest frame.
- ``photoband`` -- which band(s) to generate. In ``sinmap`` this is a
  comma-separated list (e.g. ``"0,2,5"``). In ``pipemap`` all 7 bands are
  always generated.
- ``deconv`` -- enable Richardson-Lucy deconvolution of the output map.
- ``psfdecv`` -- PSF sigma for deconvolution (arcsec). If ``0``, the PSF is
  estimated as ``sqrt(fiber_size^2 + kernel_size^2)``.


pipecube-Specific Features
--------------------------

``pipecube`` has several features beyond what ``sincube`` offers:

**Spectral chunking for large mosaics**

When a target has more than 100 exposures, the full spectral range cannot fit
in memory at once. ``pipecube`` automatically splits the cube into spectral
chunks:

- 101--500 exposures: 7 chunks (1000 A each)
- 501--1000 exposures: 13 chunks (500 A each)
- >1000 exposures: 25 chunks (250 A each)

The chunks are then co-added into a single cube when ``mergecube: True``.

**Spatial tiling**

Set ``nsp`` to a value > 1 to split the merged cube into an ``nsp x nsp``
spatial grid. This is useful for very large mosaics that exceed memory limits.

**RSS extraction**

When ``cube2rss: True``, ``pipecube`` extracts an RSS from each generated
cube and writes it to ``out_pathrss``.

**Resume support**

If a run is interrupted during spectral chunk processing, set
``index_reset`` to the chunk index where it stopped to resume without
reprocessing earlier chunks. Set ``force_merge: True`` to skip cube generation
entirely and only run the co-addition step.


pipemap-Specific Features
--------------------------

``pipemap`` generates all seven photometric bands for each target and then
automatically creates composite color JPEG images:

- **OHS** -- [OIII] (blue) + H-alpha (green) + [SII] (red)
- **OOH** -- [OII] (blue) + [OIII] (green) + H-alpha (red)
- **gri** -- g (blue) + r (green) + i (red)
- **griOOH** -- broad-band + narrow-band composite with configurable weights

If ``deconv`` is enabled, deconvolved versions of each image and composite
are also produced.


Parameter Comparison
--------------------

The table below summarizes which parameters are available in each command.

.. list-table::
   :header-rows: 1
   :widths: 20 12 12 12 12

   * - Parameter
     - sincube
     - pipecube
     - sinmap
     - pipemap
   * - ``survey_type``
     - Y
     - Y
     - Y
     - Y
   * - ``out_path``
     - Y
     - Y [#auto]_
     - Y
     - Y [#auto]_
   * - ``redux_ver``
     - Y
     - Y
     - Y
     - Y
   * - ``redux_dir``
     - Y
     - Y
     - Y
     - Y
   * - ``sigm_s``
     - Y
     -
     - Y
     -
   * - ``pix_s``
     - Y
     -
     - Y
     -
   * - ``alph_s``
     - Y
     -
     - Y
     -
   * - ``type``
     -
     - Y
     -
     - Y
   * - ``flu16``
     - Y
     - Y
     -
     -
   * - ``nameF``
     - Y
     -
     - Y
     -
   * - ``nameL``
     -
     - Y
     -
     - Y
   * - ``coord_cen``
     - Y
     -
     - Y
     -
   * - ``coord_cenL``
     -
     - Y
     -
     - Y
   * - ``use_slitmap``
     - Y
     - Y
     - Y
     - Y
   * - ``pbars``
     - Y
     - Y
     - Y
     - Y
   * - ``spec_range``
     - Y
     -
     -
     -
   * - ``path_sas``
     - Y
     - Y
     - Y
     - Y
   * - ``basename``
     - Y
     - Y
     - Y
     - Y
   * - ``errors``
     - Y
     - Y
     -
     -
   * - ``fac_sizeY``
     - Y
     - Y
     - Y
     - Y
   * - ``fac_sizeX``
     - Y
     - Y
     - Y
     - Y
   * - ``cent``
     - Y
     - Y
     - Y
     - Y
   * - ``pathF``
     - Y
     - Y
     - Y
     - Y
   * - ``mergecube``
     -
     - Y
     -
     -
   * - ``nsp``
     -
     - Y
     -
     -
   * - ``cube2rss``
     -
     - Y
     -
     -
   * - ``out_pathrss``
     -
     - Y [#auto]_
     -
     -
   * - ``nproc``
     -
     - Y
     -
     - Y
   * - ``force_merge``
     -
     - Y
     -
     -
   * - ``index_reset``
     -
     - Y
     -
     -
   * - ``redshift``
     -
     -
     - Y
     - Y
   * - ``photoband``
     -
     -
     - Y
     -
   * - ``deconv``
     -
     -
     - Y
     - Y
   * - ``psfdecv``
     -
     -
     - Y
     - Y

.. [#auto] Path has ``<redux_ver>/<type>/`` appended automatically.


Quick-Start Examples
--------------------

**Generate a single cube with custom kernel:**

.. code-block:: bash

   3dcubegen sincube -c config_sincube.yml

**Batch-generate cubes for multiple targets at type-c resolution:**

.. code-block:: bash

   3dcubegen pipecube -c config_pipecube.yml

**Generate an H-alpha narrow-band map for one target:**

.. code-block:: bash

   3dcubegen sinmap -c config_sinmap.yml

**Batch-generate all maps and color composites:**

.. code-block:: bash

   3dcubegen pipemap -c config_pipemap.yml

All commands can also be run entirely from CLI flags without a config file.
See ``3dcubegen <command> --help`` for the full list of options.
