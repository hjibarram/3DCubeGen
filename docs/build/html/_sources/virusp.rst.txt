VirusP Pipeline
===============

Overview
--------

VirusP (VIRUS-P) is a fiber-fed IFU spectrograph on the 2.7 m Harlan J. Smith
Telescope at McDonald Observatory (HET heritage instrument). It features 246
fibers with 4.16" diameter, covering a 1.7' x 1.7' field of view, operating in
the optical range (3400--6800 A typical).


CLI Entry Point
---------------

The VirusP pipeline is accessed via the ``vpcubetools`` command:

.. code-block:: bash

   vpcubetools vpcube [options]

The ``vpcube`` command generates a reconstructed data cube from one or more
VirusP observation pointings.


Basic Usage
-----------

.. code-block:: bash

   # With a config file:
   vpcubetools vpcube -c config_vpcube.yml

   # With CLI options:
   vpcubetools vpcube --name_list Name1 Name2 --redux_dir /path/to/rss \
       --out_path out_cubes/ --namef OutputName


CLI Options
-----------

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Option
     - Default
     - Description
   * - ``-c, --config_file``
     - ``''``
     - YAML configuration file (overrides CLI options)
   * - ``-i, --name_list``
     - --
     - List of input observation pointing names
   * - ``-o, --out_path``
     - ``out_cubes/``
     - Output directory for the cube
   * - ``-r, --redux_dir``
     - ``''``
     - Path to the directory containing reduced RSS files
   * - ``-s, --sigm_s``
     - ``32.2``
     - Smoothing kernel size (arcsec)
   * - ``-p, --pix_s``
     - ``0.0``
     - Spaxel size (0 = ``0.75 * sigm_s``)
   * - ``-a, --alph_s``
     - ``2.0``
     - Kernel shape factor
   * - ``-f, --flu16``
     - ``True``
     - Use 10\ :sup:`-16` erg/s/cm\ :sup:`2`/A units
   * - ``-n, --namef``
     - --
     - Output root name for the data cube
   * - ``-m, --fib``
     - ``4.16``
     - Fiber diameter (arcsec)
   * - ``-t, --base_namewcs``
     - ``NAME_wcs.txt``
     - WCS text file basename template
   * - ``-y, --basename``
     - ``NAME_abscal_HST.fits``
     - Input RSS file basename template
   * - ``-h, --basenamec``
     - ``vpCube-NAME.fits``
     - Output cube basename template
   * - ``-e, --errors``
     - ``True``
     - Enable error propagation
   * - ``-u, --fac_sizey``
     - ``1.1``
     - FoV extent factor in Y
   * - ``-x, --fac_sizex``
     - ``1.1``
     - FoV extent factor in X
   * - ``-k, --radvel``
     - ``True``
     - Apply heliocentric velocity correction
   * - ``-z, --adrcor``
     - ``True``
     - Apply Atmospheric Differential Refraction correction
   * - ``-g, --spec_range``
     - ``(None, None)``
     - Spectral range ``(lambda_min, lambda_max)`` in Angstroms
   * - ``-b, --pbars``
     - ``True``
     - Show progress bars


Configuration File
------------------

The VirusP pipeline accepts a YAML configuration file. A template is provided
at ``CubeGen/virustools/configfiles/config_vpcube.yml``:

.. code-block:: yaml

   ---
   vpcube:
     - survey_type: virusP
       name_list: [Name1, Name2]
       out_path: out_cubes/
       redux_dir: ''
       fib: 4.16
       sigm_s: 4.16
       pix_s: 0.0
       alph_s: 2.0
       flu16: True
       nameF: NameF
       spec_range: [0, 0]
       pbars: True
       base_nameWCS: NAME_wcs.txt
       basename: NAME_abscal_HST.fits
       basenameC: vpCube-NAME.fits
       errors: True
       fac_sizeY: 1.1
       fac_sizeX: 1.1
       radvel: True
       adrcor: True


VirusP-Specific Features
-------------------------

**Separate WCS files**

Unlike LVM (which reads astrometry from the SLITMAP extension of each RSS
file), VirusP uses external WCS text files. The ``base_nameWCS`` parameter
specifies the basename template; for each pointing in ``name_list``, the code
reads ``<redux_dir>/<name>_wcs.txt`` to obtain the fiber positions on sky.

**Heliocentric velocity correction**

When ``radvel: True``, the pipeline applies a heliocentric velocity correction
to shift the wavelength solution from the topocentric to the heliocentric
frame.

**ADR correction**

When ``adrcor: True``, fiber positions are adjusted per wavelength channel to
correct for Atmospheric Differential Refraction. This activates the ADR-aware
interpolation mode described in :doc:`algorithm`, where the kernel weights
become wavelength-dependent.

**Configurable fiber diameter**

The ``fib`` parameter sets the fiber diameter (default 4.16"), which affects
both the radius cutoff and the flux area scaling factor. This differs from LVM
(fixed at 35.3") and MEGARA (derived from the bundle geometry).


Differences from LVM
---------------------

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Aspect
     - VirusP
     - LVM
   * - Fiber count
     - 246
     - ~1800 science fibers
   * - Fiber diameter
     - 4.16" (configurable)
     - 35.3" (fixed)
   * - Default ``sigm_s``
     - 4.16" (fiber-matched)
     - 32.2"
   * - Astrometry source
     - External WCS text files
     - SLITMAP extension in RSS FITS
   * - Directory structure
     - Flat (``redux_dir/``)
     - SAS tree (version/tilegroup/tile/mjd)
   * - Heliocentric correction
     - Built-in (``radvel``)
     - Applied during RSS reduction
   * - ADR correction
     - Built-in (``adrcor``)
     - Not available
