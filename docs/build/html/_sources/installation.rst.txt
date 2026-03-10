Installation & Getting Started
==============================

Prerequisites
-------------

- Python 3.11 or later (3.13 recommended)
- A conda or virtual environment is recommended to avoid dependency conflicts

Installing from Source
----------------------

Clone the repository and install in development mode:

.. code-block:: bash

   git clone -b dev https://github.com/hjibarram/3DCubeGen.git
   cd 3DCubeGen
   pip install -e .

Verify the installation by checking that the three CLI entry points are
available:

.. code-block:: bash

   3dcubegen --help      # LVM pipeline
   pipemegara --help     # MEGARA pipeline
   vpcubetools --help    # VirusP tools


Core Dependencies
-----------------

The following packages are installed automatically by ``pip install``:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Package
     - Purpose
   * - ``astropy``
     - FITS I/O, WCS transformations, coordinate handling
   * - ``numpy``
     - Array operations throughout the codebase
   * - ``scipy``
     - Spatial distance calculations, interpolation
   * - ``matplotlib``
     - Plotting utilities
   * - ``pyyaml``
     - YAML configuration file parsing
   * - ``click`` / ``cloup``
     - Command-line interface framework
   * - ``tqdm``
     - Progress bars during cube reconstruction
   * - ``sdeconv``
     - Richardson-Lucy deconvolution for map sharpening


Optional Dependencies
---------------------

These are not installed automatically and are only needed for specific
pipelines:

**MEGARA pipeline** -- requires the ``numina`` framework and ``megaradrp``
(MEGARA Data Reduction Pipeline):

.. code-block:: bash

   pip install numina megaradrp

**Deconvolution** -- the ``scikit-image`` package provides additional
deconvolution routines:

.. code-block:: bash

   pip install scikit-image


Quick Verification
------------------

After installation, a quick way to confirm everything works is to generate the
help text for each command:

.. code-block:: bash

   # LVM commands
   3dcubegen sincube --help
   3dcubegen pipecube --help
   3dcubegen sinmap --help
   3dcubegen pipemap --help

   # MEGARA commands
   pipemegara rundrp --help
   pipemegara runfulldrp --help

   # VirusP commands
   vpcubetools vpcube --help
