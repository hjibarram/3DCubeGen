YAML Configuration Reference
============================

3DCubeGen uses YAML configuration files to control LVM cube generation.
Two command modes are available: ``sincube`` for generating a single cube, and
``pipecube`` for batch-processing multiple targets through a pipeline.

Template configuration files are provided in the ``Configfiles/`` directory:

- ``Configfiles/config_sincube.yml`` -- single-cube mode
- ``Configfiles/config_pipecube.yml`` -- pipeline mode


sincube Configuration
---------------------

The ``sincube`` section defines parameters for generating a single data cube
from one target.

.. list-table::
   :header-rows: 1
   :widths: 18 10 12 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``survey_type``
     - string
     - ``LVM``
     - Survey identifier. Must contain ``LVM`` for LVM data.
   * - ``out_path``
     - string
     - ``out_cubes/``
     - Output directory for the generated cube FITS files.
   * - ``redux_ver``
     - string
     - ``1.1.1.dev0``
     - DRP reduction version of the input ``lvmCFrame`` files.
   * - ``redux_dir``
     - string
     - ``''``
     - Path to the redux directory. When empty, defaults to
       ``<path_sas>/sdsswork/lvm/spectro/redux``.
   * - ``sigm_s``
     - float
     - ``32.2``
     - Smoothing kernel size (arcsec) used to reconstruct the cube.
       This sets the spatial resolution.
   * - ``pix_s``
     - float
     - ``0.0``
     - Output spaxel size (arcsec). When set to ``0``, it is auto-calculated
       as ``0.75 * sigm_s``.
   * - ``alph_s``
     - float
     - ``2.0``
     - Kernel shape factor. A value of ``2.0`` produces a Gaussian kernel.
       Values larger than 2 produce a sharper (more top-hat-like) kernel.
   * - ``flu16``
     - bool
     - ``True``
     - If ``True``, output cube flux units are 10\ :sup:`-16` erg/s/cm\ :sup:`2`/A.
       If ``False``, units are erg/s/cm\ :sup:`2`/A.
   * - ``nameF``
     - string
     - --
     - Root name for the target. Used as both the output cube filename stem
       and the exposure list filename (read from ``pathF + nameF``).
   * - ``use_slitmap``
     - bool
     - ``True``
     - Use the astrometry stored in the slitmap header of the RSS files.
   * - ``pbars``
     - bool
     - ``True``
     - Enable or disable progress bars during processing.
   * - ``spec_range``
     - list
     - ``[0, 0]``
     - Spectral range ``[lambda_min, lambda_max]`` in Angstroms for the
       output cube. Set to ``[0, 0]`` to use the full spectral range from the
       input RSS files.
   * - ``path_sas``
     - string
     - ``''``
     - Path to the local SDSS SAS directory. Used only if the
       ``SAS_BASE_DIR`` environment variable is not defined (see
       :ref:`sas-env-var`).
   * - ``basename``
     - string
     - ``lvmSFrame-NAME.fits``
     - Template filename of the input RSS files. ``NAME`` is replaced at
       runtime with the exposure number.
   * - ``errors``
     - bool
     - ``True``
     - Enable error propagation. When ``True``, the output cube includes an
       error extension computed via inverse-variance weighting.
   * - ``fac_sizeY``
     - float
     - ``1.1``
     - Multiplicative factor applied to the field-of-view extent in the Y
       (declination) direction.
   * - ``fac_sizeX``
     - float
     - ``1.1``
     - Multiplicative factor applied to the field-of-view extent in the X
       (right ascension) direction.
   * - ``cent``
     - bool
     - ``True``
     - Enable centering of the output field-of-view on the coordinates given
       by ``coord_cen``.
   * - ``coord_cen``
     - list
     - ``[0, 0]``
     - Central coordinates ``[RA, Dec]`` in degrees for the output cube.
       Used when ``cent`` is ``True``.
   * - ``pathF``
     - string
     - ``''``
     - Directory path where the exposure list file is located. The code opens
       ``pathF + nameF``, so this must end with ``/`` if non-empty.


sincube Example
^^^^^^^^^^^^^^^

.. code-block:: yaml

   ---
   sincube:
     - survey_type: LVM
       out_path: out_cubes/
       redux_ver: 1.1.1.dev0
       redux_dir: ''
       sigm_s: 32.2
       pix_s: 0.0
       alph_s: 2
       flu16: True
       nameF: NGC_1234
       use_slitmap: True
       pbars: True
       spec_range: [0, 0]
       path_sas: ''
       basename: lvmSFrame-NAME.fits
       errors: True
       fac_sizeY: 1.1
       fac_sizeX: 1.1
       cent: True
       coord_cen: [180.0, -30.0]
       pathF: exposure_lists/


pipecube Configuration
----------------------

The ``pipecube`` section defines parameters for batch-processing multiple
targets. It shares most parameters with ``sincube`` but adds pipeline-specific
options and uses a predefined kernel type instead of explicit ``sigm_s``/``pix_s``
values.

.. list-table::
   :header-rows: 1
   :widths: 18 10 12 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``survey_type``
     - string
     - ``LVM``
     - Survey identifier. Must contain ``LVM``.
   * - ``out_path``
     - string
     - ``out_cubes/``
     - Base output directory. The pipeline automatically appends
       ``<redux_ver>/<type>/`` to this path and creates the directory if it
       does not exist.
   * - ``redux_ver``
     - string
     - ``1.1.1.dev0``
     - DRP reduction version.
   * - ``redux_dir``
     - string
     - ``''``
     - Path to the redux directory. When empty, defaults to
       ``<path_sas>/sdsswork/lvm/spectro/redux``.
   * - ``type``
     - string
     - ``c``
     - Kernel type code that sets the spatial resolution. See
       :ref:`kernel-types` below.
   * - ``flu16``
     - bool
     - ``True``
     - Flux unit scaling (same as sincube).
   * - ``nameL``
     - list
     - --
     - List of target names, e.g. ``[NGC_1234, NGC_5678]``. Each name is used
       to find the corresponding exposure list file at ``pathF + name``.
   * - ``use_slitmap``
     - bool
     - ``True``
     - Use slitmap astrometry.
   * - ``pbars``
     - bool
     - ``True``
     - Enable progress bars.
   * - ``path_sas``
     - string
     - ``''``
     - Local SAS directory path (see :ref:`sas-env-var`).
   * - ``basename``
     - string
     - ``lvmSFrame-NAME.fits``
     - RSS file template name.
   * - ``errors``
     - bool
     - ``True``
     - Enable error propagation.
   * - ``fac_sizeY``
     - float
     - ``1.1``
     - FoV Y extent factor.
   * - ``fac_sizeX``
     - float
     - ``1.1``
     - FoV X extent factor.
   * - ``cent``
     - bool
     - ``True``
     - Enable FoV centering.
   * - ``coord_cenL``
     - list
     - ``[[0,0]]``
     - List of ``[RA, Dec]`` coordinate pairs (degrees), one per target in
       ``nameL``. Used when ``cent`` is ``True``.
   * - ``pathF``
     - string
     - ``''``
     - Directory containing the exposure list files.
   * - ``mergecube``
     - bool
     - ``False``
     - Enable cube merging. When a target has more than 100 exposures, the
       cube is built in spectral chunks and then co-added.
   * - ``nsp``
     - int
     - ``0``
     - Number of spatial tiles per axis for the merged cube (``nsp x nsp``
       grid). Set to ``0`` or ``1`` to disable spatial tiling.
   * - ``cube2rss``
     - bool
     - ``True``
     - After cube generation, extract an RSS from the cube.
   * - ``out_pathrss``
     - string
     - ``out_rss/``
     - Base output directory for extracted RSS files. The pipeline appends
       ``<redux_ver>/<type>/`` automatically.
   * - ``nproc``
     - int
     - ``3``
     - Number of parallel threads for the kernel interpolation.
   * - ``force_merge``
     - bool
     - ``False``
     - When ``True``, skip cube generation and only run the co-addition step.
       Useful for restarting after a partial run.
   * - ``index_reset``
     - int
     - ``0``
     - Starting index for spectral chunk processing. Set to a non-zero value
       to resume from a specific chunk after interruption.


pipecube Example
^^^^^^^^^^^^^^^^

.. code-block:: yaml

   ---
   pipecube:
     - survey_type: LVM
       out_path: out_cubes/
       redux_ver: 1.1.1.dev0
       redux_dir: ''
       type: c
       flu16: True
       nameL: [NGC_1234, NGC_5678]
       use_slitmap: True
       pbars: True
       path_sas: ''
       basename: lvmSFrame-NAME.fits
       errors: True
       fac_sizeY: 1.1
       fac_sizeX: 1.1
       cent: True
       coord_cenL: [[180.0, -30.0], [200.0, -25.0]]
       pathF: exposure_lists/
       mergecube: False
       nsp: 0
       cube2rss: True
       out_pathrss: out_rss/
       nproc: 3
       force_merge: False
       index_reset: 0


.. _kernel-types:

Kernel Type Codes
-----------------

The ``pipecube`` command uses a single-letter ``type`` code to select a
predefined kernel size. The kernel size (``sigm_s``) and spaxel size
(``pix_s``) are computed as:

.. code-block:: none

   sigm_s = 8.8 * 32 * multiplier    (arcsec)
   pix_s  = 8.8 * 0.75 * 32 * multiplier  (arcsec)

where the multiplier depends on the type code:

.. list-table::
   :header-rows: 1
   :widths: 8 12 20 20

   * - Code
     - Multiplier
     - sigm_s (arcsec)
     - pix_s (arcsec)
   * - ``a``
     - 16
     - 4505.6
     - 3379.2
   * - ``b``
     - 8
     - 2252.8
     - 1689.6
   * - ``c``
     - 4
     - 1126.4
     - 844.8
   * - ``d``
     - 2
     - 563.2
     - 422.4
   * - ``e``
     - 1
     - 281.6
     - 211.2
   * - ``f``
     - 1/2
     - 140.8
     - 105.6
   * - ``g``
     - 1/4
     - 70.4
     - 52.8
   * - ``h``
     - 1/8
     - 35.2
     - 26.4
   * - ``i``
     - 1/16
     - 17.6
     - 13.2
   * - ``j``
     - 1/32
     - 8.8
     - 6.6

Larger multipliers produce coarser spatial resolution (bigger kernels and
spaxels), suitable for large mosaics. Smaller multipliers give finer resolution
for individual pointings.

.. note::

   For targets with fewer than 3 exposures, types ``i`` and ``j`` are
   automatically upgraded to type ``h`` to avoid undersampling.


.. _sas-env-var:

SAS Directory Resolution
-------------------------

The code resolves the SAS base directory in this order:

1. If the ``SAS_BASE_DIR`` environment variable is set, its value is used
   (with a trailing ``/`` appended).
2. Otherwise, the ``path_sas`` config parameter is used (with a trailing
   ``/`` appended).

The redux directory is then constructed as:

.. code-block:: none

   <SAS_BASE_DIR>/sdsswork/lvm/spectro/redux

unless ``redux_dir`` is explicitly set to a non-empty value.


Output Path Construction
------------------------

**sincube**: Output files are written directly to ``out_path``.

**pipecube**: The pipeline appends ``<redux_ver>/<type>/`` to both
``out_path`` and ``out_pathrss``. For example, with:

.. code-block:: yaml

   out_path: out_cubes/
   out_pathrss: out_rss/
   redux_ver: 1.1.1.dev0
   type: c

the effective output paths become:

.. code-block:: none

   out_cubes/1.1.1.dev0/c/
   out_rss/1.1.1.dev0/c/

These directories are created automatically if they do not exist.
