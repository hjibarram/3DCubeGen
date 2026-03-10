Exposure List File Format
=========================

The exposure list is a simple CSV file that tells 3DCubeGen which RSS
(Row-Stacked Spectra) frames to combine into a data cube. It is read by the
``read_explist()`` function in ``CubeGen/tools/tools.py``.


File Format
-----------

The file is a headerless, comma-separated text file with **3 columns**:

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Column
     - Type
     - Description
   * - 1
     - integer
     - ``tile_id`` -- the SDSS-V tile identifier for this pointing
   * - 2
     - integer
     - ``mjd`` -- Modified Julian Date of the observation
   * - 3
     - string
     - ``exposure_number`` -- the exposure number (e.g. ``00012345``)

**Parsing rules:**

- Lines containing ``#`` are treated as comments and skipped.
- Lines containing ``TILE`` are skipped (allows an optional header line).
- Leading and trailing whitespace is stripped from each field.
- Empty fields after splitting are removed.


Example File
------------

.. code-block:: text

   # Exposure list for NGC_1234
   # TILE, MJD, EXPN
   10842,60310,00003456
   10842,60310,00003457
   10842,60310,00003458
   10842,60311,00003510
   10843,60311,00003511


Tile Group Derivation
---------------------

Internally, ``read_explist()`` derives a **tile group** string from each
``tile_id``. This tile group is used to locate the RSS files in the SDSS
directory tree:

- For most tiles: the first 4 characters of the tile ID followed by ``XX``.
  For example, tile ``10842`` becomes group ``1084XX``.
- Special case: tile ``11111`` maps to group ``0011XX``.


How the Exposure List Connects to Config Parameters
----------------------------------------------------

The exposure list filename and location are controlled by config parameters:

**sincube mode:**

- ``nameF`` -- the filename of the exposure list (also used as the output
  cube name stem)
- ``pathF`` -- the directory containing the file

The code opens ``pathF + nameF``, so if ``pathF`` is non-empty it must end
with ``/``.

**pipecube mode:**

- ``nameL`` -- a list of filenames (one per target)
- ``pathF`` -- the directory containing all the files

For each target ``nameL[i]``, the code opens ``pathF + nameL[i]``.


RSS File Path Resolution
------------------------

Each row in the exposure list maps to an RSS file on disk. The full path is
constructed from the config parameters and the parsed exposure list fields:

.. code-block:: none

   <redux_dir>/<redux_ver>/<tile_group>/<tile_id>/<mjd>/<basename>

where ``<basename>`` has ``NAME`` replaced with the exposure number. For
example, with these settings:

.. code-block:: yaml

   redux_dir: /data/sdsswork/lvm/spectro/redux
   redux_ver: 1.1.1.dev0
   basename: lvmSFrame-NAME.fits

and an exposure list row:

.. code-block:: text

   10842,60310,00003456

the resolved RSS file path would be:

.. code-block:: none

   /data/sdsswork/lvm/spectro/redux/1.1.1.dev0/1084XX/10842/60310/lvmSFrame-00003456.fits
