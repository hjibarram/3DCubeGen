Kernel Smoothing Parameters
===========================

The ``sincube``/``sinmap`` and ``pipecube``/``pipemap`` commands use the same
underlying kernel interpolation algorithm but specify its smoothing parameters
in different ways.  This page explains both conventions and shows how to
translate between them.


sincube / sinmap: explicit parameters
--------------------------------------

The single-target commands accept three kernel parameters directly in the
config file:

.. list-table::
   :header-rows: 1
   :widths: 12 55 33

   * - Parameter
     - Description
     - Default
   * - ``sigm_s``
     - Smoothing kernel size :math:`\sigma` in arcsec.  Controls the spatial
       resolution of the output cube: a larger value smooths over more fibers
       and produces coarser effective resolution.
     - ``32.2``
   * - ``pix_s``
     - Output spaxel size in arcsec.  Set to ``0`` to auto-compute as
       ``0.75 * sigm_s``.
     - ``0.0``
   * - ``alph_s``
     - Kernel shape factor :math:`\alpha`.  ``2.0`` gives a Gaussian profile;
       values above 2.0 produce a sharper, more top-hat-like kernel.
     - ``2.0``

Example ``sincube`` config section::

    sincube:
      - sigm_s: 35.2
        pix_s: 0.0      # auto → 0.75 × 35.2 = 26.4 arcsec
        alph_s: 2.0


pipecube / pipemap: type letter
--------------------------------

The pipeline commands replace the three parameters above with a single
``type`` letter (``a`` through ``j``).  The letter is passed to
``tools.kernel_pipe()`` which computes ``sigm_s`` and ``pix_s`` from a fixed
formula based on the LVM fiber size:

.. math::

   \sigma   &= \frac{17.6}{2} \times 32 \times v = 8.8 \times 32 \times v \\[4pt]
   p        &= 0.75\,\sigma

where :math:`v` is a scale factor that doubles with each step from ``j`` to
``a``.  The shape factor :math:`\alpha` is **always 2.0** (Gaussian) in
pipeline mode and cannot be changed through the config.

Example ``pipecube`` config section::

    pipecube:
      - type: h


.. _kernel-types:

Type-to-parameter mapping
--------------------------

The table below gives the ``sigm_s`` and ``pix_s`` values produced by each
``type`` letter, together with the equivalent ``sincube`` config block.

.. list-table::
   :header-rows: 1
   :widths: 7 8 18 18 10 39

   * - ``type``
     - scale :math:`v`
     - ``sigm_s`` (arcsec)
     - ``pix_s`` (arcsec)
     - ``alph_s``
     - Equivalent ``sincube`` / ``sinmap`` values
   * - ``a``
     - 16
     - 4505.6
     - 3379.2
     - 2.0
     - ``sigm_s: 4505.6``, ``pix_s: 3379.2``, ``alph_s: 2.0``
   * - ``b``
     - 8
     - 2252.8
     - 1689.6
     - 2.0
     - ``sigm_s: 2252.8``, ``pix_s: 1689.6``, ``alph_s: 2.0``
   * - ``c``
     - 4
     - 1126.4
     - 844.8
     - 2.0
     - ``sigm_s: 1126.4``, ``pix_s: 844.8``, ``alph_s: 2.0``
   * - ``d``
     - 2
     - 563.2
     - 422.4
     - 2.0
     - ``sigm_s: 563.2``, ``pix_s: 422.4``, ``alph_s: 2.0``
   * - ``e``
     - 1
     - 281.6
     - 211.2
     - 2.0
     - ``sigm_s: 281.6``, ``pix_s: 211.2``, ``alph_s: 2.0``
   * - ``f``
     - 1/2
     - 140.8
     - 105.6
     - 2.0
     - ``sigm_s: 140.8``, ``pix_s: 105.6``, ``alph_s: 2.0``
   * - ``g``
     - 1/4
     - 70.4
     - 52.8
     - 2.0
     - ``sigm_s: 70.4``, ``pix_s: 52.8``, ``alph_s: 2.0``
   * - ``h``
     - 1/8
     - 35.2
     - 26.4
     - 2.0
     - ``sigm_s: 35.2``, ``pix_s: 26.4``, ``alph_s: 2.0``
   * - ``i``
     - 1/16
     - 17.6
     - 13.2
     - 2.0
     - ``sigm_s: 17.6``, ``pix_s: 13.2``, ``alph_s: 2.0``
   * - ``j``
     - 1/32
     - 8.8
     - 6.6
     - 2.0
     - ``sigm_s: 8.8``, ``pix_s: 6.6``, ``alph_s: 2.0``

The scale factor :math:`v` doubles for each step up the alphabet, so the
resolution coarsens by a factor of two per step.  The ``sigm_s`` base value of
281.6 arcsec (type ``e``, :math:`v = 1`) equals ``8.8 × 32``, where 8.8 arcsec
is the LVM fiber radius.

.. note::

   When ``pix_s: 0`` is used in ``sincube``/``sinmap``, the auto-computed
   spaxel size is ``0.75 × sigm_s`` -- the same ratio used by ``kernel_pipe()``.
   Setting ``pix_s: 0`` in ``sincube`` is therefore equivalent to the pipeline
   behaviour.


Differences beyond the kernel
------------------------------

Even with identical ``sigm_s``, ``pix_s``, and ``alph_s`` values, the two
modes differ in a few other respects:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Feature
     - ``sincube`` / ``sinmap``
     - ``pipecube`` / ``pipemap``
   * - ``alph_s`` configurable
     - Yes
     - No (hardcoded to 2.0)
   * - Targets per run
     - One (``nameF``)
     - Many (``nameL``)
   * - RSS extraction
     - Not available
     - Optional (``cube2rss: True``)
   * - Spectral chunking
     - Not available
     - Automatic for > 100 exposures
   * - Cube co-addition
     - Not available
     - Optional (``mergecube: True``)
   * - Output path
     - ``out_path`` as given
     - ``out_path/<redux_ver>/<type>/``


Small-exposure guard in pipecube
---------------------------------

When fewer than three exposures are available for a target, ``pipecube``
automatically promotes types ``i`` and ``j`` up to type ``h``:

.. code-block:: text

   < 3 exposures  +  type i or j  →  forced to type h
   (sigm_s = 35.2 arcsec,  pix_s = 26.4 arcsec)

This prevents the finest-resolution kernels from being used when there is
insufficient fiber coverage to support them.  There is no equivalent guard in
``sincube``; if you set ``sigm_s`` to 17.6 or 8.8 arcsec with very few
exposures the kernel will be applied as specified.
