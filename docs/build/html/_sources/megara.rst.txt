MEGARA Pipeline
===============

Overview
--------

`MEGARA <https://www.gtc.iac.es/instruments/megara/>`_ (Multi-Espectrografo en
GTC de Alta Resolucion para Astronomia) is an optical IFU and multi-object
spectrograph on the 10.4 m Gran Telescopio Canarias (GTC). It features a
hexagonal fiber bundle of 623 fibers covering a 12.5" x 11.3" field of view,
with fiber core diameter of 0.62".

3DCubeGen provides a complete pipeline for MEGARA that integrates the
``numina`` / ``megaradrp`` Data Reduction Pipeline (DRP) with the kernel
interpolation cube reconstruction.


CLI Entry Point
---------------

The MEGARA pipeline is accessed via the ``pipemegara`` command:

.. code-block:: bash

   pipemegara <command> [options]

Two commands are available:

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Command
     - Purpose
   * - ``rundrp``
     - Manual step-by-step DRP execution + cube reconstruction
   * - ``runfulldrp``
     - Fully automated pipeline: sort observations, reduce calibrations,
       calibrate science data, and reconstruct cubes


rundrp
------

The ``rundrp`` command runs the MEGARA DRP reduction steps individually and
then reconstructs a data cube. It requires the ``numina`` framework to be
installed and the MEGARA DRP configuration files (``obsresult-*.yaml``,
``requirements-*.yaml``) to be present in the working directory.

The reduction follows this calibration chain:

1. **Bias** -- master bias creation
2. **Traces** -- fiber trace identification
3. **Wavelength calibration** -- arc lamp wavelength solution
4. **Fiber flat** -- fiber-to-fiber response correction
5. **Standard star** -- spectrophotometric calibration
6. **Science** -- final science frame reduction + cube reconstruction

.. code-block:: bash

   pipemegara rundrp --redux_dir /path/to/data --stdar_t HD12345 --vph LR-R

Key options:

.. list-table::
   :header-rows: 1
   :widths: 20 12 68

   * - Option
     - Default
     - Description
   * - ``-o, --out_path``
     - ``out_cubes``
     - Output directory for the final cube
   * - ``-r, --redux_dir``
     - ``''``
     - Path to the working directory with raw data
   * - ``-s, --sigm_s``
     - ``0.35``
     - Smoothing kernel size (arcsec)
   * - ``-p, --pix_s``
     - ``0.0``
     - Spaxel size (0 = same as ``sigm_s``)
   * - ``-a, --alph_s``
     - ``2.0``
     - Kernel shape factor
   * - ``-m, --stdar_t``
     - --
     - Name of the spectrophotometric standard star
   * - ``-y, --vph``
     - ``LR-R``
     - VPH grating name (e.g., LR-B, LR-R, LR-Z, MR-R, HR-R)
   * - ``-u, --aper_std``
     - ``4.0``
     - Aperture radius for spectrophotometric calibration (arcsec)
   * - ``-b, --basename``
     - ``megCube-NAME.fits``
     - Output cube basename template


runfulldrp
----------

The ``runfulldrp`` command automates the entire reduction workflow. It:

1. Sorts observation files by type (bias, arc, flat, standard, science)
2. Generates the DRP configuration files automatically
3. Runs all calibration steps
4. Performs spectrophotometric calibration iteratively
5. Reconstructs the final data cube
6. Extracts a 1D spectrum

It can be driven by CLI options or by a YAML configuration file.

.. code-block:: bash

   pipemegara runfulldrp -c config_runfulldrp.yml

   # Or with CLI options:
   pipemegara runfulldrp --run_name GTC5-24AMEX --run_number 0001 \
       --obrun_path /path/to/raw --redux_path /path/to/redux \
       --out_path /path/to/output

Key additional options:

.. list-table::
   :header-rows: 1
   :widths: 22 14 64

   * - Option
     - Default
     - Description
   * - ``-t, --run_name``
     - ``GTC5-24AMEX``
     - Name of the GTC observation run
   * - ``-g, --run_number``
     - ``0001``
     - Observation block number (string format)
   * - ``-m, --obrun_path``
     - ``''``
     - Path to the raw observation data
   * - ``-r, --redux_path``
     - ``''``
     - Working directory for intermediate reduction products
   * - ``-c, --config_file``
     - ``''``
     - Path to YAML config file (overrides CLI options)
   * - ``-e, --forcear``
     - ``False``
     - Force use of ThAr arc lamp for wavelength calibration


Configuration File
------------------

The ``runfulldrp`` command accepts a YAML configuration file. A template is
provided at ``CubeGen/megaratools/configfiles/config_runfulldrp.yml``:

.. code-block:: yaml

   ---
   run_drp:
     - pipe_type: FULL
       run_name: GTC5-24AMEX
       run_number: '0001'
       out_path: out_cubes
       redux_ver: 1.1.1.dev0
       redux_path: ''
       obrun_path: ''
       sigm_s: 0.35
       pix_s: 0.0
       alph_s: 2.0
       flu16: True
       name: None
       basename: megCube-NAME.fits
       fac_sizey: 1.0
       fac_sizex: 1.0
       aper_std: 4.0
       forcear: False


Differences from LVM
---------------------

While MEGARA uses the same core kernel interpolation algorithm, several
defaults and behaviors differ:

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Aspect
     - MEGARA
     - LVM
   * - Default ``sigm_s``
     - 0.35" (sub-arcsecond)
     - 32.2" (wide fibers)
   * - ``pix_s`` when 0
     - Same as ``sigm_s``
     - ``0.75 * sigm_s``
   * - Fiber geometry
     - 623 fibers, hexagonal, 0.62" core
     - ~1800 science fibers, 35.3" pitch
   * - VPH selection
     - Multiple gratings (LR-B/R/Z, MR, HR)
     - N/A (fixed spectral range)
   * - Spectrophotometric cal
     - Built into pipeline (iterative)
     - External (pre-calibrated SFrames)
   * - DRP integration
     - Runs ``numina`` + ``megaradrp``
     - Reads pre-reduced RSS frames
