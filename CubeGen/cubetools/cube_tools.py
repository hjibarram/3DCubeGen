from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import CubeGen.megaratools.megtools as mtools
import CubeGen.tools.tools as tools
import numpy as np
from tqdm.notebook import tqdm
from tqdm import tqdm as tqdmT
from PIL import Image


def astromatch(file0,file1,sig=2,plotview=False):
    """
    Compare the astrometric registration of two reconstructed data cubes.

    The function collapses each input cube along the spectral axis to
    generate a 2D integrated-flux image. A two-dimensional PSF model is
    then fitted to each image using :func:`evaluate_2dPSF` in
    ``CubeGen.megaratools.megtools``.

    The fitted centroid of the reference cube (``file0``) is converted
    to sky coordinates using its WCS and subsequently projected onto the
    pixel coordinate system of ``file1``. This allows the measured
    centroid displacement between the two cubes to be compared with the
    displacement expected from their WCS solutions.

    Parameters
    ----------
    file0 : str or path-like
        Path to the reference FITS data cube. The primary HDU must
        contain a three-dimensional array with dimensions
        ``(wavelength, x, y)`` and a valid celestial WCS.

    file1 : str or path-like
        Path to the FITS data cube whose astrometric registration is
        compared with ``file0``. The primary HDU must contain a
        three-dimensional array and a valid celestial WCS.

    sig : float, optional
        Initial or characteristic sigma used by
        :func:`CubeGen.megaratools.megtools.evaluate_2dPSF`.
        The default is 2.

    plotview : bool, optional
        If True, display the integrated images and fitted PSF models.

    Returns
    -------
    dx : float
        Difference in the measured PSF centroid along the x pixel
        coordinate,

        ``x0 - x1``.

    dy : float
        Difference in the measured PSF centroid along the y pixel
        coordinate,

        ``y0 - y1``.

    dx_wcs : float
        Difference between the measured x centroid in ``file1`` and
        the expected position of the ``file0`` centroid after projection
        through the WCS of ``file1``.

    dy_wcs : float
        Difference between the measured y centroid in ``file1`` and
        the expected position of the ``file0`` centroid after projection
        through the WCS of ``file1``.

    Notes
    -----
    The integrated images are generated as

    .. math::

        I(x,y) = \\sum_\\lambda F(\\lambda,x,y),

    ignoring NaN values.

    Two types of offsets are calculated:

    1. The direct centroid displacement between the two cubes,

       .. math::

           \\Delta x = x_0 - x_1,

           \\Delta y = y_0 - y_1.

    2. The displacement relative to the WCS prediction. The centroid
       measured in ``file0`` is transformed to sky coordinates and then
       projected onto the pixel grid of ``file1``,

       .. math::

           \\Delta x_{\\mathrm{WCS}} = x_1 - x_{0\\rightarrow1},

           \\Delta y_{\\mathrm{WCS}} = y_1 - y_{0\\rightarrow1}.

    A value close to zero for ``dx_wcs`` and ``dy_wcs`` indicates that
    the relative position of the source is consistent with the WCS
    solutions of the two cubes.

    Examples
    --------
    >>> from CubeGen.cubetools.cube_tools import astromatch
    >>> dx, dy, dx_wcs, dy_wcs = astromatch(
    ...     "cube_reference.fits",
    ...     "cube_target.fits",
    ...     sig=2
    ... )
    >>> print(dx_wcs, dy_wcs)

    See Also
    --------
    CubeGen.megaratools.megtools.evaluate_2dPSF
        Fits a two-dimensional PSF model and determines its centroid.
    """
    [spec0, hdr0]=fits.getdata(file0, 0, header=True)
    
    [spec1, hdr1]=fits.getdata(file1, 0, header=True)

    # Collapse the cubes along the spectral axis.
    try:
        nz0,nx0,ny0=spec0.shape
        map0=np.nansum(spec0,axis=0)
    except:
        nx0,ny0=spec0.shape
        map0=np.copy(spec0)
    try:
        nz1,nx1,ny1=spec1.shape
        map1=np.nansum(spec1,axis=0)
    except:
        nx1,ny1=spec1.shape
        map1=np.copy(spec1)
    print(file0)
    print(file1)

    # Determine the PSF centroid in each reconstructed image.
    x0,y0,ds_m0,psf0,model0=mtools.evaluate_2dPSF(map0,model=True,sig=sig,plotview=plotview)
    x1,y1,ds_m1,psf1,model1=mtools.evaluate_2dPSF(map1,model=True,sig=sig,plotview=plotview)

    print("x_0=",x0,"y_0=",y0,"sigma_0=",ds_m0,"psf_0=",psf0)
    print("x_1=",x1,"y_1=",y1,"sigma_1=",ds_m1,"psf_1=",psf1)
    
    # Extract the celestial component of the WCS.
    wcs0 = WCS(hdr0)
    wcs0=wcs0.celestial
    wcs1 = WCS(hdr1)
    wcs1=wcs1.celestial
    
    # Convert measured centroids to celestial coordinates.
    sky0=pixel_to_skycoord(x0,y0,wcs0)
    sky1=pixel_to_skycoord(x1,y1,wcs1)
    val0=sky0.to_string('hmsdms')
    val1=sky1.to_string('hmsdms')
    print(val0,'RA0,DEC0')
    print(val1,'RA1,DEC1')
    
    # Expected location of the source from file0 in the pixel system
    # of file1.
    xpos,ypos=skycoord_to_pixel(sky0,wcs1)
    print("x_1_sk=",xpos,"y_1_sk=",ypos)
    
    dx = x0 - x1
    dy = y0 - y1

    dx_wcs = x1 - xpos
    dy_wcs = y1 - ypos

    print("Dx=", dx, dx_wcs)
    print("Dy=", dy, dy_wcs)
    return dx, dy, dx_wcs, dy_wcs


def crop_image(names, cube, dir1='./', dir2='./', dir3='./', apt='_gri'):
    """
    Reproject imaging data onto the spatial grid of an IFU data cube.

    This function takes a set of FITS images and resamples them onto the
    spatial pixel grid defined by a reference three-dimensional data cube.
    The transformation between the imaging data and the IFU cube is
    performed using their celestial World Coordinate System (WCS)
    solutions.

    Before reprojection, the image flux is corrected for the difference
    in pixel area between the input image and the reference cube. The
    reprojected images are stored as a multi-plane FITS file.

    When three input images are supplied, the function also constructs
    an RGB representation using an astronomical magnitude scaling and
    saves it as a JPEG image.

    Parameters
    ----------
    names : sequence of str
        Base names of the input FITS images, without the ``.fits``
        extension. Each file is expected to contain a two-dimensional
        image with a valid celestial WCS.

        For RGB generation, three images are expected. In the current
        implementation they are interpreted in the order supplied and
        mapped directly onto the three RGB channels.

    cube : str
        Filename of the reference IFU cube. The primary HDU must contain
        a three-dimensional data array and a valid celestial WCS.

    dir1 : str, optional
        Directory containing the input imaging FITS files.
        The default is ``'.'``.

    dir2 : str, optional
        Directory containing the reference IFU cube.
        The default is ``'.'``.

    apt : str, optional
        Suffix added to the output filenames. The default is ``'_gri'``.

    Returns
    -------
    None
        The function writes its results directly to disk.

    Outputs
    -------
    FITS file
        A multi-plane FITS image containing the reprojected input images.
        The spatial WCS is inherited from the reference IFU cube.

        The output filename is constructed from ``cube`` by replacing
        ``'.fits.gz'`` with ``apt`` and appending ``'.fits'``.

    JPEG file
        An RGB image generated from the first three reprojected images.
        The output filename follows the same convention as the FITS
        product, with the ``'.jpeg'`` extension.

    Notes
    -----
    The spatial pixel scale of the reference cube is estimated from the
    CD matrix as

    .. math::

        s_x = 3600
        \\sqrt{\\mathrm{CD1\\_1}^2 + \\mathrm{CD1\\_2}^2}

    and

    .. math::

        s_y = 3600
        \\sqrt{\\mathrm{CD2\\_1}^2 + \\mathrm{CD2\\_2}^2},

    where the resulting scales are expressed in arcseconds per pixel.

    The corresponding pixel area is

    .. math::

        A = s_x s_y.

    Input images are rescaled according to the ratio between the
    reference-cube pixel area, :math:`A_0`, and the input-image pixel
    area, :math:`A_1`,

    .. math::

        I' = I \\frac{A_0}{A_1}.

    Each output pixel is then transformed from the reference-cube pixel
    grid to celestial coordinates and subsequently into the pixel
    coordinate system of the input image. The image intensity at the
    resulting position is evaluated using ``map_interpolB``.

    For the RGB output, the image values are converted to an
    astronomical magnitude-like scale using

    .. math::

        m = -2.5 \\log_{10}
        \\left(
        I \\frac{3.631\\times10^{-6}}{f_0}
        \\right),

    where ``f0`` is defined by the ``zerop`` array in the function.

    The RGB generation currently assumes exactly three input images.

    Examples
    --------
    >>> from CubeGen.cubetools.cube_tools import crop_image
    >>>
    >>> images = [
    ...     "frame-i-003366-5-0074",
    ...     "frame-r-003366-5-0074",
    ...     "frame-g-003366-5-0074",
    ... ]
    >>>
    >>> crop_image(
    ...     images,
    ...     "J1338+4816_R.fits.gz",
    ...     dir1="J1338+4816/",
    ...     dir2="/data/cubes/",
    ...     apt="_gri",
    ... )

    See Also
    --------
    astromatch
        Compare the astrometric registration of two reconstructed cubes.
    """

    file1 = dir2 + cube

    # Read reference IFU cube.
    spec1, hdr1 = fits.getdata(file1, 0, header=True)

    # Determine the reference pixel scale and pixel area.
    dx = np.sqrt(hdr1['CD1_1']**2.0 +
                 hdr1['CD1_2']**2.0) * 3600.0
    dy = np.sqrt(hdr1['CD2_1']**2.0 +
                 hdr1['CD2_2']**2.0) * 3600.0
    A0 = dx * dy

    

    # Integrated image of the reference cube.
    try:
        nz1, nx1, ny1 = spec1.shape
        map1 = np.nansum(spec1, axis=0)
    except:
        nx1, ny1 = spec1.shape
        map1= np.copy(spec1)

    # Celestial WCS of the IFU cube.
    wcs1 = WCS(hdr1).celestial

    nt = len(names)

    # Reprojected imaging cube.
    pdl_cube_new = np.zeros([nt, nx1, ny1])

    ct = 0

    for name in names:

        fig1 = np.zeros([nx1, ny1])

        # Read external imaging data.
        cube_file = dir1 + name + '.fits'
        pdl_cube, hdr0 = fits.getdata(
            cube_file, 0, header=True
        )

        # Determine input-image pixel area.
        dx = np.sqrt(hdr0['CD1_1']**2.0 +
                     hdr0['CD1_2']**2.0) * 3600.0
        dy = np.sqrt(hdr0['CD2_1']**2.0 +
                     hdr0['CD2_2']**2.0) * 3600.0

        A1 = dx * dy

        # Celestial WCS of input image.
        wcs0 = WCS(hdr0).celestial

        # Correct for differences in pixel area.
        pdl_cube = pdl_cube * A0 / A1

        # Reproject every IFU spatial pixel onto the input image.
        for i in range(0, nx1):
            for j in range(0, ny1):

                sky1 = pixel_to_skycoord(j, i, wcs1)

                xpos, ypos = skycoord_to_pixel(
                    sky1, wcs0
                )

                val = tools.map_interpolB(
                    pdl_cube,
                    ypos,
                    xpos
                )

                fig1[i, j] = val

        pdl_cube_new[ct, :, :] = fig1
        ct += 1

    # ---------------------------------------------------------
    # Write reprojected imaging cube
    # ---------------------------------------------------------

    h1 = fits.PrimaryHDU(pdl_cube_new)
    h_k = h1.header

    keys = list(hdr1.keys())

    for key in keys:
        h_k[key] = hdr1[key]
        h_k.comments[key] = hdr1.comments[key]

    try:
        # Remove the original spectral-axis WCS.
        del h_k['CDELT3']
        del h_k['CRPIX3']
        del h_k['CRVAL3']
    except:
        print('2D map')

    h_k['BUNIT'] = hdr0['BUNIT']

    hlist = fits.HDUList([h1])
    hlist.writeto(
        dir3 + cube.replace('.fits.gz', apt) + '.fits',
        overwrite=True
    )

    # ---------------------------------------------------------
    # Generate RGB representation
    # ---------------------------------------------------------

    zerop = [3730.0, 3730.0, 3631.0]

    rgb_cube = np.zeros(
        [nx1, ny1, 3],
        dtype="uint8"
    )

    pdl_img_r = -2.5 * np.log10(
        pdl_cube_new[1, :, :] *
        3.631e-6 / zerop[1]
    )

    vmax = (
        np.amin(
            pdl_img_r[
                np.where(pdl_img_r > 0)
            ]
        ) - 0.5
    )

    vmin = 25.0

    for i in range(0, 3):

        pdl_img = -2.5 * np.log10(
            pdl_cube_new[i, :, :] *
            3.631e-6 / zerop[i]
        )

        pdl_img = (
            (pdl_img - vmin) /
            (vmax - vmin) *
            256
        )

        pdl_img[pdl_img < 0] = 0
        pdl_img[pdl_img > 255] = 255

        rgb_cube[:, :, i] = np.flipud(pdl_img)


    im = Image.fromarray(rgb_cube)

    im.save(
        dir3 + cube.replace('.fits.gz', apt) + '.jpeg',
        quality=100
    )



def coad_cube(name, dir1='', dir2='', vphs=None, patch=True, verbose=False,
              pbars=True, notebook=True):
    """
    Co-add reconstructed IFU datacubes from multiple spectral bands into
    a single wavelength-continuous datacube.

    The function reads reconstructed datacubes from ``dir1`` using
    ``name`` as the common filename and the spectral-band or VPH
    identifiers specified in ``vphs``. The individual cubes are spatially
    registered using their WCS information and resampled onto a common
    wavelength grid.

    Fluxes in overlapping wavelength regions are combined and, when
    required, multiplicative corrections are applied to match adjacent
    spectral bands.

    The resulting co-added cube is written to ``dir2`` as
    ``name.fits.gz``. When ``patch=True``, an additional
    ``name_patch.fits.gz`` file containing the correction-factor map
    is generated.

    Parameters
    ----------
    name : str
        Common base filename of the input datacubes, excluding the VPH
        suffix and ``.fits.gz`` extension.

        For example, if::

            name = 'MK883'
            vphs = ['B', 'G', 'R']

        the function searches for::

            MK883_B.fits.gz
            MK883_G.fits.gz
            MK883_R.fits.gz

    dir1 : str, optional
        Directory containing the input reconstructed datacubes.
        Default is ``''``, corresponding to the current directory.

    dir2 : str, optional
        Directory where the final co-added datacube and optional patch
        map are written. Default is ``''``, corresponding to the current
        directory.

    vphs : list of str or None, optional
        Ordered list containing the identifiers of the spectral bands
        or VPHs to be co-added. At least two bands must be provided and
        a maximum of three bands is currently supported.

        The order of the entries must follow increasing wavelength.
        For a three-band configuration, the first, second, and third
        entries correspond to the blue, intermediate, and red spectral
        ranges, respectively.

        For example::

            vphs = ['B', 'G', 'R']

        or::

            vphs = ['LR-B', 'LR-V', 'LR-R']

        The identifiers are used to construct the input filenames as::

            <name>_<vph>.fits.gz

        If ``None``, the function prints a message requesting the VPH
        definition and returns without processing the cubes.

    patch : bool, optional
        If True, calculate and save the multiplicative correction-factor
        map used to match overlapping spectral bands. Default is True.

    verbose : bool, optional
        If True, print additional information about the processing steps

    pbar : bool, optional
        If True, display a progress bar during the spatial loop. Default is True.

    notebook : bool, optional
        If True, use the Jupyter notebook version of the progress bar.

    Returns
    -------
    None
        The function writes the resulting FITS products to disk. If the
        VPH list is not defined or fewer than two valid spectral bands
        are available, the function returns without generating an output
        cube.

    Notes
    -----
    The input datacubes are expected to contain a three-dimensional
    primary HDU with wavelength along the third FITS axis.

    The wavelength solution is obtained from ``CRPIX3``, ``CRVAL3`` and
    either ``CD3_3`` or ``CDELT3`` in the FITS headers.

    The spectral bands supplied in ``vphs`` must be ordered from shorter
    to longer wavelengths.

    The output FITS file contains:

    * Primary HDU
        Co-added flux datacube.
    * ``Error_cube``
        Estimated uncertainty datacube.
    * ``BADPIXELMASK``
        Integer bad-pixel mask.
    * ``match_F_cube``
        Spatial map of the multiplicative matching factor.

    Examples
    --------
    Co-add three spectral bands::

        >>> coad_cube(
        ...     'MK883',
        ...     dir1='/data/cubes/',
        ...     dir2='/data/coadd/',
        ...     vphs=['B', 'G', 'R']
        ... )

    This searches for::

        /data/cubes/MK883_B.fits.gz
        /data/cubes/MK883_G.fits.gz
        /data/cubes/MK883_R.fits.gz

    and generates::

        /data/coadd/MK883.fits.gz
        /data/coadd/MK883_patch.fits.gz

    Different VPH identifiers can be specified, for example::

        >>> coad_cube(
        ...     'MK883',
        ...     dir1='/data/cubes/',
        ...     dir2='/data/coadd/',
        ...     vphs=['LR-B', 'LR-V', 'LR-R']
        ... )

    If the VPH list is omitted::

        >>> coad_cube('MK883')
        Define at least two VPHs/bands.
    """
    # ------------------------------------------------------------------
    # Validate VPH/band definition
    # ------------------------------------------------------------------
    if vphs is None:
        print("Define at least two VPHs/bands.")
        return

    if len(vphs) < 2:
        print("At least two VPHs/bands must be provided.")
        return

    if len(vphs) > 3:
        print("A maximum of three VPHs/bands is currently supported.")
        return

    # Ensure directories end with "/".
    if dir1 != '' and not dir1.endswith('/'):
        dir1 += '/'

    if dir2 != '' and not dir2.endswith('/'):
        dir2 += '/'

    # Band identifiers.
    vph_B = vphs[0]

    if len(vphs) >= 2:
        vph_G = vphs[1]
    else:
        vph_G = None

    if len(vphs) >= 3:
        vph_R = vphs[2]
    else:
        vph_R = None

    # ------------------------------------------------------------------
    # Read available spectral-band cubes
    # ------------------------------------------------------------------

    try:
        cube_file = dir1 + name + '_' + vph_B + '.fits.gz'
        specB, hdrB = fits.getdata(cube_file, 0, header=True)

        nz0, nx0, ny0 = specB.shape
        band_B = True

    except:
        band_B = False

    if vph_G is not None:

        try:
            cube_file = dir1 + name + '_' + vph_G + '.fits.gz'
            specG, hdrG = fits.getdata(cube_file, 0, header=True)

            ## Original empirical flux correction.
            #specG = specG * 1.833

            nz1, nx1, ny1 = specG.shape
            band_G = True

        except:
            band_G = False

    else:
        band_G = False

    if vph_R is not None:

        try:
            cube_file = dir1 + name + '_' + vph_R + '.fits.gz'
            specR, hdrR = fits.getdata(cube_file, 0, header=True)

            nz2, nx2, ny2 = specR.shape
            band_R = True

        except:
            band_R = False

    else:
        band_R = False

    # ------------------------------------------------------------------
    # Determine first and last available spectral bands
    # ------------------------------------------------------------------

    if band_B:
        nx = nx0
        ny = ny0
        init = 1

    elif band_G:
        nx = nx1
        ny = ny1
        init = 2

    elif band_R:
        nx = nx2
        ny = ny2
        init = 3

    else:
        print("No valid input cubes were found.")
        return

    if band_R:
        fint = 3

    elif band_G:
        fint = 2

    elif band_B:
        fint = 1

    else:
        print("No valid input cubes were found.")
        return

    # ------------------------------------------------------------------
    # Construct wavelength arrays and celestial WCS
    # ------------------------------------------------------------------

    if band_B:

        wcsB = WCS(hdrB).celestial

        crpix = hdrB["CRPIX3"]

        try:
            cdelt = hdrB["CD3_3"]
        except:
            cdelt = hdrB["CDELT3"]

        crval = hdrB["CRVAL3"]

        waveB = crval + cdelt * (
            np.arange(nz0) + 1 - crpix
        )

        pixsB = cdelt

    if band_G:

        wcsG = WCS(hdrG).celestial

        crpix = hdrG["CRPIX3"]

        try:
            cdelt = hdrG["CD3_3"]
        except:
            cdelt = hdrG["CDELT3"]

        crval = hdrG["CRVAL3"]

        waveG = crval + cdelt * (
            np.arange(nz1) + 1 - crpix
        )

        pixsG = cdelt

    if band_R:

        wcsR = WCS(hdrR).celestial

        crpix = hdrR["CRPIX3"]

        try:
            cdelt = hdrR["CD3_3"]
        except:
            cdelt = hdrR["CDELT3"]

        crval = hdrR["CRVAL3"]

        waveR = crval + cdelt * (
            np.arange(nz2) + 1 - crpix
        )

        pixsR = cdelt

    # ------------------------------------------------------------------
    # Select reference cube
    # ------------------------------------------------------------------

    if init == 1:

        wcs = wcsB
        min_wave = np.nanmin(waveB)

        hdr0 = hdrB

        nxi = nx0
        nyi = ny0
        nzi = nz0

    if init == 2:

        wcs = wcsG
        min_wave = np.nanmin(waveG)

        hdr0 = hdrG

        nxi = nx1
        nyi = ny1
        nzi = nz1

    if init == 3:
        print("Only one spectral band is available. Exiting.")
        return

    if fint == 3:

        max_wave = np.nanmax(waveR)

    if fint == 2:

        if init == 2:
            print("Only one spectral band is available. Exiting.")
            return

        max_wave = np.nanmax(waveG)

    if fint == 1:
        print("Only one spectral band is available. Exiting.")
        return

    # ------------------------------------------------------------------
    # Determine common wavelength sampling
    # ------------------------------------------------------------------

    if band_B and band_R:

        cdelt = np.amax(
            np.array([pixsB, pixsR])
        )

        pix_mB = int(np.ceil(cdelt / pixsB))
        pix_mR = int(np.ceil(cdelt / pixsR))

    if band_B and band_G:

        cdelt = np.amax(
            np.array([pixsB, pixsG])
        )

        pix_mB = int(np.round(cdelt / pixsB))
        pix_mG = int(np.round(cdelt / pixsG))

    if band_G and band_R:

        cdelt = np.amax(
            np.array([pixsG, pixsR])
        )

        pix_mG = int(np.round(cdelt / pixsG))
        pix_mR = int(np.round(cdelt / pixsR))

    if band_B and band_G and band_R:

        cdelt = np.amax(
            np.array([pixsB, pixsG, pixsR])
        )

        pix_mB = int(np.round(cdelt / pixsB))
        pix_mG = int(np.round(cdelt / pixsG))
        pix_mR = int(np.round(cdelt / pixsR))

    n_pix = int(
        np.round(
            (max_wave - min_wave) / cdelt
        )
    )

    crval = min_wave
    crpix = 1

    waveF = crval + cdelt * (
        np.arange(n_pix) + 1 - crpix
    )

    # ------------------------------------------------------------------
    # Allocate output cubes
    # ------------------------------------------------------------------

    IFU_coadd = np.zeros(
        [n_pix, nx, ny]
    )

    IFU_coaddE = np.zeros(
        [n_pix, nx, ny]
    )

    IFU_coaddB = np.zeros(
        [n_pix, nx, ny],
        dtype=int
    )

    patch_map = np.zeros(
        [nx, ny]
    )

    # ------------------------------------------------------------------
    # Spatial loop
    # ------------------------------------------------------------------
    if pbars:
        if notebook:
            pbar=tqdm(total=nx)
        else:     
            pbar=tqdmT(total=nx)  
    for i in range(nx):
        for j in range(ny):

            temp_spec = np.copy(
                IFU_coadd[:, i, j]
            )

            pix1 = i
            pix2 = j

            sky1 = pixel_to_skycoord(
                pix2,
                pix1,
                wcs
            )

            val1 = sky1.to_string(
                'hmsdms'
            )
            if verbose:
                print(
                    val1,
                    'RA,DEC'
                )

            ypos0, xpos0 = skycoord_to_pixel(
                sky1,
                wcs
            )

            xpos0 = int(
                np.round(xpos0)
            )

            ypos0 = int(
                np.round(ypos0)
            )

            # ----------------------------------------------------------
            # Reference spectral band
            # ----------------------------------------------------------

            if (
                xpos0 >= 0 and
                xpos0 <= nxi - 1 and
                ypos0 >= 0 and
                ypos0 <= nyi - 1
            ):

                if init == 2:

                    spec1 = specG[
                        :, xpos0, ypos0
                    ]

                    if pix_mG > 1:
                        spec1 = tools.median_a(
                            spec1,
                            lw=pix_mG
                        )

                else:

                    spec0 = specB[
                        :, xpos0, ypos0
                    ]

                    if pix_mB > 1:
                        spec0 = tools.median_a(
                            spec0,
                            lw=pix_mB * 2
                        )
                if verbose:
                    print(
                        xpos0,
                        ypos0,
                        'POS0',
                        i,
                        j
                    )

            else:

                if init == 2:
                    spec1 = np.zeros(nzi)

                else:
                    spec0 = np.zeros(nzi)

            # ----------------------------------------------------------
            # Intermediate spectral band
            # ----------------------------------------------------------

            if band_G and init == 1:

                ypos1a, xpos1a = skycoord_to_pixel(
                    sky1,
                    wcsG
                )

                xpos1 = int(
                    np.round(xpos1a)
                )

                ypos1 = int(
                    np.round(ypos1a)
                )

                if (
                    xpos1 >= 0 and
                    xpos1 <= nx1 - 1 and
                    ypos1 >= 0 and
                    ypos1 <= ny1 - 1
                ):

                    spec1 = specG[
                        :, xpos1, ypos1
                    ]

                    if np.nansum(spec1) != 0:

                        spec1 = tools.cube_interpolB(
                            specG,
                            xpos1a,
                            ypos1a
                        )

                    if pix_mG > 1:

                        spec1 = tools.median_a(
                            spec1,
                            lw=pix_mG
                        )
                    if verbose:
                        print(
                            xpos1,
                            ypos1,
                            'POS1',
                            i,
                            j
                        )

                else:

                    spec1 = np.zeros(nz1)

            # ----------------------------------------------------------
            # Red spectral band
            # ----------------------------------------------------------

            if band_R:

                ypos2a, xpos2a = skycoord_to_pixel(
                    sky1,
                    wcsR
                )

                # Original empirical spatial offsets.
                #ypos2a = ypos2a + 0.3999999999999986
                #xpos2a = xpos2a - 0.1999999999999993 - 0.5

                xpos2 = int(
                    np.round(xpos2a)
                )

                ypos2 = int(
                    np.round(ypos2a)
                )

                if (
                    xpos2 >= 0 and
                    xpos2 <= nx2 - 1 and
                    ypos2 >= 0 and
                    ypos2 <= ny2 - 1
                ):

                    spec2 = specR[
                        :, xpos2, ypos2
                    ]

                    if np.nansum(spec2) != 0:

                        spec2 = tools.cube_interpolB(
                            specR,
                            xpos2a,
                            ypos2a
                        )

                    if pix_mR > 1:

                        spec2 = tools.median_a(
                            spec2,
                            lw=pix_mR
                        )
                    if verbose:
                        print(
                            xpos2a,
                            ypos2a,
                            'POS2',
                            i,
                            j
                        )
                        print(
                            xpos2, ypos2,
                            'POS2_ROUND',
                            i, j
                        )

                else:

                    spec2 = np.zeros(nz2)

            # ----------------------------------------------------------
            # B + R
            # ----------------------------------------------------------

            if band_B and band_R and not band_G:

                nt1 = np.where(
                    waveB >= np.nanmin(waveR)
                )

                nt2 = np.where(
                    waveR <= np.nanmax(waveB)
                )

                nt1i = np.where(
                    waveB <= np.nanmin(waveR) + 13
                )

                nt2s = np.where(
                    waveR >= np.nanmax(waveB) - 13
                )

                ntF = np.where(
                    (waveF <= np.nanmax(waveB)) &
                    (waveF >= np.nanmin(waveR))
                )

                ntFs = np.where(
                    waveF >= np.nanmax(waveB) - 13
                )

                ntFi = np.where(
                    waveF <= np.nanmin(waveR) + 13
                )

                fc = 1.0

                if len(nt1[0]) > 0:

                    specBF = np.interp(
                        waveF[ntF],
                        waveB[nt1],
                        spec0[nt1],
                        left=0.,
                        right=0.
                    )

                    specRF = np.interp(
                        waveF[ntF],
                        waveR[nt2],
                        spec2[nt2],
                        left=0.,
                        right=0.
                    )

                    Bflux = np.nanmean(specBF)
                    Rflux = np.nanmean(specRF)

                    ft = Rflux / Bflux

                    if 0.1 <= ft <= 5.5:
                        fc = ft
                    else:
                        fc = 1.0

                    if patch:
                        patch_map[i, j] = fc
                    if verbose:
                        print(
                            "factor=",
                            fc
                        )

                    specBRF = (
                        specBF * fc + specRF
                    ) / 2.0

                    temp_spec[ntF] = specBRF

                if len(nt1i[0]) > 0:

                    specBFi = np.interp(
                        waveF[ntFi],
                        waveB[nt1i],
                        spec0[nt1i],
                        left=0.,
                        right=0.
                    )

                    temp_spec[ntFi] = (
                        specBFi * fc
                    )

                if len(nt2s[0]) > 0:

                    specRFs = np.interp(
                        waveF[ntFs],
                        waveR[nt2s],
                        spec2[nt2s],
                        left=0.,
                        right=0.
                    )

                    temp_spec[ntFs] = specRFs

            # ----------------------------------------------------------
            # B + G
            # ----------------------------------------------------------

            if band_B and band_G and not band_R:

                nt1 = np.where(
                    waveB >= np.nanmin(waveG)
                )

                nt2 = np.where(
                    waveG <= np.nanmax(waveB)
                )

                nt1i = np.where(
                    waveB <= np.nanmin(waveG) + 3
                )

                nt2s = np.where(
                    waveG >= np.nanmax(waveB) - 3
                )

                ntF = np.where(
                    (waveF <= np.nanmax(waveB)) &
                    (waveF >= np.nanmin(waveG))
                )

                ntFs = np.where(
                    waveF >= np.nanmax(waveB) - 3
                )

                ntFi = np.where(
                    waveF <= np.nanmin(waveG) + 3
                )

                fc = 1.0

                if len(nt1[0]) > 0:

                    specBF = np.interp(
                        waveF[ntF],
                        waveB[nt1],
                        spec0[nt1],
                        left=0.,
                        right=0.
                    )

                    specGF = np.interp(
                        waveF[ntF],
                        waveG[nt2],
                        spec1[nt2],
                        left=0.,
                        right=0.
                    )

                    Bflux = np.nanmean(specBF)
                    Gflux = np.nanmean(specGF)

                    ft = Gflux / Bflux

                    if 0.1 <= ft <= 15.5:
                        fc = ft
                    else:
                        fc = 1.0

                    if patch:
                        patch_map[i, j] = fc
                    if verbose:
                        print(
                            "factor=",
                            fc
                        )

                    specBGF = (
                        specBF * fc + specGF
                    ) / 2.0

                    temp_spec[ntF] = specBGF

                if len(nt1i[0]) > 0:

                    specBFi = np.interp(
                        waveF[ntFi],
                        waveB[nt1i],
                        spec0[nt1i],
                        left=0.,
                        right=0.
                    )

                    temp_spec[ntFi] = (
                        specBFi * fc
                    )

                if len(nt2s[0]) > 0:

                    specGFs = np.interp(
                        waveF[ntFs],
                        waveG[nt2s],
                        spec1[nt2s],
                        left=0.,
                        right=0.
                    )

                    temp_spec[ntFs] = specGFs

            # ----------------------------------------------------------
            # G + R
            # ----------------------------------------------------------

            if band_G and band_R and not band_B:

                nt1 = np.where(
                    waveG >= np.nanmin(waveR)
                )

                nt2 = np.where(
                    waveR <= np.nanmax(waveG)
                )

                nt1i = np.where(
                    waveG <= np.nanmin(waveR) + 1
                )

                nt2s = np.where(
                    waveR >= np.nanmax(waveG) - 1
                )

                ntF = np.where(
                    (waveF <= np.nanmax(waveG)) &
                    (waveF >= np.nanmin(waveR))
                )

                ntFs = np.where(
                    waveF >= np.nanmax(waveG) - 1
                )

                ntFi = np.where(
                    waveF <= np.nanmin(waveR) + 1
                )

                if len(nt1[0]) > 0:

                    specGF = np.interp(
                        waveF[ntF],
                        waveG[nt1],
                        spec1[nt1],
                        left=0.,
                        right=0.
                    )

                    specRF = np.interp(
                        waveF[ntF],
                        waveR[nt2],
                        spec2[nt2],
                        left=0.,
                        right=0.
                    )

                    specGRF = (
                        specGF + specRF
                    ) / 2.0

                    temp_spec[ntF] = specGRF

                if len(nt1i[0]) > 0:

                    specGFi = np.interp(
                        waveF[ntFi],
                        waveG[nt1i],
                        spec1[nt1i],
                        left=0.,
                        right=0.
                    )

                    temp_spec[ntFi] = specGFi

                if len(nt2s[0]) > 0:

                    specRFs = np.interp(
                        waveF[ntFs],
                        waveR[nt2s],
                        spec2[nt2s],
                        left=0.,
                        right=0.
                    )

                    temp_spec[ntFs] = specRFs

            # ----------------------------------------------------------
            # B + G + R
            # ----------------------------------------------------------

            if band_B and band_G and band_R:

                nt1a = np.where(
                    waveB >= np.nanmin(waveG)
                )

                nt2a = np.where(
                    waveG <= np.nanmax(waveB)
                )

                ntFa = np.where(
                    (waveF <= np.nanmax(waveB)) &
                    (waveF >= np.nanmin(waveG))
                )

                fc = 1.0

                if len(nt1a[0]) > 0:

                    specBFa = np.interp(
                        waveF[ntFa],
                        waveB[nt1a],
                        spec0[nt1a],
                        left=0.,
                        right=0.
                    )

                    specGFa = np.interp(
                        waveF[ntFa],
                        waveG[nt2a],
                        spec1[nt2a],
                        left=0.,
                        right=0.
                    )

                    Bflux = np.nanmean(specBFa)
                    Gflux = np.nanmean(specGFa)

                    ft = Bflux / Gflux

                    if 0.01 <= ft <= 5.0:
                        fc = ft
                    else:
                        fc = 1.0

                    if patch:
                        patch_map[i, j] = fc
                    if verbose:
                        print(
                            "factorGB=",
                            fc
                        )

                    specBGFa = (
                        specGFa * fc + specBFa
                    ) / 2.0

                    temp_spec[ntFa] = specBGFa

                nt1b = np.where(
                    waveB >= np.nanmin(waveR)
                )

                nt2b = np.where(
                    waveR <= np.nanmax(waveB)
                )

                ntFb = np.where(
                    (waveF <= np.nanmax(waveB)) &
                    (waveF >= np.nanmin(waveR))
                )

                if len(nt1b[0]) > 0:

                    specBFb = np.interp(
                        waveF[ntFb],
                        waveB[nt1b],
                        spec0[nt1b],
                        left=0.,
                        right=0.
                    )

                    specRFb = np.interp(
                        waveF[ntFb],
                        waveR[nt2b],
                        spec2[nt2b],
                        left=0.,
                        right=0.
                    )

                    specBRFb = (
                        specBFb + specRFb
                    ) / 2.0

                    temp_spec[ntFb] = specBRFb

                nt1c = np.where(
                    waveG >= np.nanmin(waveR)
                )

                nt2c = np.where(
                    waveR <= np.nanmax(waveG)
                )

                ntFc = np.where(
                    (waveF <= np.nanmax(waveG)) &
                    (waveF >= np.nanmin(waveR))
                )

                if len(nt1c[0]) > 0:

                    specGFc = np.interp(
                        waveF[ntFc],
                        waveG[nt1c],
                        spec1[nt1c],
                        left=0.,
                        right=0.
                    )

                    specRFc = np.interp(
                        waveF[ntFc],
                        waveR[nt2c],
                        spec2[nt2c],
                        left=0.,
                        right=0.
                    )

                    specGRFc = (
                        specGFc + specRFc
                    ) / 2.0

                    temp_spec[ntFc] = specGRFc

                nt1bc = np.where(
                    (waveB <= np.nanmin(waveG) + 1) &
                    (waveB <= np.nanmin(waveR) + 1)
                )

                ntFbc = np.where(
                    (waveF <= np.nanmin(waveG) + 1) &
                    (waveF <= np.nanmin(waveR) + 1)
                )

                if len(nt1bc[0]) > 0:

                    specB_GRFbc = np.interp(
                        waveF[ntFbc],
                        waveB[nt1bc],
                        spec0[nt1bc],
                        left=0.,
                        right=0.
                    )

                    temp_spec[ntFbc] = (
                        specB_GRFbc# * fc
                    )

                nt1ac = np.where(
                    (waveG <= np.nanmin(waveR) + 1) &
                    (waveG >= np.nanmax(waveB))
                )

                ntFac = np.where(
                    (waveF <= np.nanmin(waveR) + 1) &
                    (waveF >= np.nanmax(waveB))
                )

                if len(nt1ac[0]) > 0:

                    specG_BRFac = np.interp(
                        waveF[ntFac],
                        waveG[nt1ac],
                        spec1[nt1ac],
                        left=0.,
                        right=0.
                    )

                    temp_spec[ntFac] = specG_BRFac * fc

                nt1ab = np.where(
                    (waveR >= np.nanmax(waveB)) &
                    (waveR >= np.nanmax(waveG))
                )

                ntFab = np.where(
                    (waveF >= np.nanmax(waveB)) &
                    (waveF >= np.nanmax(waveG))
                )

                if len(nt1ab[0]) > 0:

                    specR_BGFab = np.interp(
                        waveF[ntFab],
                        waveR[nt1ab],
                        spec2[nt1ab],
                        left=0.,
                        right=0.
                    )

                    temp_spec[ntFab] = specR_BGFab

            # ----------------------------------------------------------
            # Store combined spectrum
            # ----------------------------------------------------------

            IFU_coadd[:, i, j] = temp_spec

            # ----------------------------------------------------------
            # Estimate uncertainty
            # ----------------------------------------------------------

            nt_z = np.where(
                temp_spec != 0
            )

            if len(nt_z[0]) > 0:

                temp_specM = tools.conv(
                    temp_spec,
                    ke=5
                )

                temp_specE = np.abs(
                    temp_spec - temp_specM
                )

                temp_specE = np.sqrt(
                    tools.conv(
                        temp_specE**2.0,
                        ke=50
                    )
                )

                IFU_coaddE[:, i, j] = (
                    temp_specE * 0.1
                )

                ntp = np.where(
                    temp_spec == 0
                )

                if len(ntp[0]) > 0:

                    IFU_coaddE[
                        ntp, i, j
                    ] = 0.002

            else:

                IFU_coaddE[
                    :, i, j
                ] = 1.0

            # ----------------------------------------------------------
            # Bad-pixel mask
            # ----------------------------------------------------------

            nt_z = np.where(
                temp_spec == 0
            )

            if len(nt_z[0]) > 0:

                IFU_coaddB[
                    nt_z, i, j
                ] = 1
        if pbars:
            pbar.update(1)
    if pbars:
        pbar.close()        
    # ------------------------------------------------------------------
    # Generate output FITS file
    # ------------------------------------------------------------------

    dx = 0
    dy = 0

    h1 = fits.PrimaryHDU(
        IFU_coadd
    )

    h2 = fits.ImageHDU(
        IFU_coaddE
    )

    h3 = fits.ImageHDU(
        IFU_coaddB
    )

    h4 = fits.ImageHDU(
        patch_map
    )

    keys = list(
        hdr0.keys()
    )

    # ------------------------------------------------------------------
    # Primary HDU
    # ------------------------------------------------------------------

    h_k = h1.header

    for key in keys:

        h_k[key] = hdr0[key]
        h_k.comments[key] = (
            hdr0.comments[key]
        )

    h_k['CDELT3'] = cdelt
    h_k['CRPIX3'] = crpix
    h_k['CRVAL3'] = crval

    h_k['CRPIX1'] = (
        h_k['CRPIX1'] + dx
    )

    h_k['CRPIX2'] = (
        h_k['CRPIX2'] + dy
    )

    # ------------------------------------------------------------------
    # Error cube
    # ------------------------------------------------------------------

    h_t = h2.header

    for key in keys:

        h_t[key] = hdr0[key]
        h_t.comments[key] = (
            hdr0.comments[key]
        )

    h_t['EXTNAME'] = 'Error_cube'

    h_t['CDELT3'] = cdelt
    h_t['CRPIX3'] = crpix
    h_t['CRVAL3'] = crval

    h_t['CRPIX1'] = (
        h_t['CRPIX1'] + dx
    )

    h_t['CRPIX2'] = (
        h_t['CRPIX2'] + dy
    )

    # ------------------------------------------------------------------
    # Bad-pixel mask
    # ------------------------------------------------------------------

    h_r = h3.header

    for key in keys:

        h_r[key] = hdr0[key]
        h_r.comments[key] = (
            hdr0.comments[key]
        )

    h_r['EXTNAME'] = 'BADPIXELMASK'

    h_r['CDELT3'] = cdelt
    h_r['CRPIX3'] = crpix
    h_r['CRVAL3'] = crval

    h_r['CRPIX1'] = (
        h_r['CRPIX1'] + dx
    )

    h_r['CRPIX2'] = (
        h_r['CRPIX2'] + dy
    )

    # ------------------------------------------------------------------
    # Matching-factor map
    # ------------------------------------------------------------------

    h_w = h4.header

    for key in keys:

        h_w[key] = hdr0[key]
        h_w.comments[key] = (
            hdr0.comments[key]
        )

    h_w['EXTNAME'] = 'match_F_cube'

    for key in (
        'CDELT3',
        'CRPIX3',
        'CRVAL3'
    ):

        if key in h_w:
            del h_w[key]

    h_w['CRPIX1'] = (
        h_w['CRPIX1'] + dx
    )

    h_w['CRPIX2'] = (
        h_w['CRPIX2'] + dy
    )

    # ------------------------------------------------------------------
    # Write combined cube
    # ------------------------------------------------------------------

    output_file = (
        dir2 + name + '.fits'
    )

    hlist = fits.HDUList(
        [h1, h2, h3, h4]
    )

    hlist.update_extend()

    hlist.writeto(
        output_file,
        overwrite=True
    )

    tools.sycall(
        'gzip -f ' + output_file
    )

    # ------------------------------------------------------------------
    # Write standalone patch map
    # ------------------------------------------------------------------

    if patch:

        hp = fits.PrimaryHDU(
            patch_map
        )

        h_k = hp.header

        for key in keys:

            h_k[key] = hdr0[key]

            h_k.comments[key] = (
                hdr0.comments[key]
            )

        for key in (
            'CDELT3',
            'CRPIX3',
            'CRVAL3'
        ):

            if key in h_k:
                del h_k[key]

        h_k['CRPIX1'] = (
            h_k['CRPIX1'] + dx
        )

        h_k['CRPIX2'] = (
            h_k['CRPIX2'] + dy
        )

        patch_file = (
            dir2 +
            name +
            '_patch.fits'
        )

        hlist = fits.HDUList(
            [hp]
        )

        hlist.update_extend()

        hlist.writeto(
            patch_file,
            overwrite=True
        )

        tools.sycall(
            'gzip -f ' + patch_file
        )


def crop_cube(file0, file1, file2, spsample_copy=False,
              fac_sizeX=1.0, fac_sizeY=1.0, dx=0, dy=0,
              pbars=True, notebook=True):
    """
    Crop a datacube using the field of view of a reference datacube.

    The spatial field of view (FoV) and astrometry of ``file0`` are used
    to define the celestial region extracted from ``file1``. The output
    can either preserve the spatial sampling of ``file1`` or be
    resampled onto the spatial grid defined by ``file0``.

    The FoV can be enlarged independently along the two spatial axes
    using ``fac_sizeX`` and ``fac_sizeY``. In addition, the reference
    astrometry can be shifted by ``dx`` and ``dy`` pixels before defining
    the output FoV.

    Parameters
    ----------
    file0 : str or path-like
        Reference FITS datacube. Its celestial WCS, spatial sampling,
        orientation, and FoV define the region to extract from ``file1``.

    file1 : str or path-like
        FITS datacube to crop. The primary HDU must contain the flux cube
        with dimensions ``(wavelength, x, y)``. Extension 1 is assumed
        to contain the corresponding uncertainty cube.

    file2 : str or path-like
        Output filename. The ``.fits`` or ``.fits.gz`` extension may be
        included or omitted.

    spsample_copy : bool, optional
        If False, preserve the native spatial sampling, orientation, and
        WCS of ``file1`` and perform a simple spatial crop.

        If True, resample ``file1`` onto the spatial grid defined by
        ``file0`` using ``tools.map_interpolB``. In this mode the output
        inherits the spatial pixel scale and orientation of the reference
        cube. Default is False.

    fac_sizeX : float, optional
        Multiplicative factor applied to the FoV along the NumPy second
        spatial dimension (FITS axis 2). A value of 1.0 reproduces the
        reference FoV, while values greater than 1 enlarge it around its
        centre. Default is 1.0.

    fac_sizeY : float, optional
        Multiplicative factor applied to the FoV along the NumPy third
        spatial dimension (FITS axis 1). A value of 1.0 reproduces the
        reference FoV, while values greater than 1 enlarge it around its
        centre. Default is 1.0.

    dx : float, optional
        Astrometric displacement along the FITS x direction (axis 1),
        expressed in pixels of the reference cube. Positive values shift
        the reference FoV toward increasing x pixel coordinates.
        Default is 0.

    dy : float, optional
        Astrometric displacement along the FITS y direction (axis 2),
        expressed in pixels of the reference cube. Positive values shift
        the reference FoV toward increasing y pixel coordinates.
        Default is 0.

    pbar : bool, optional
        If True, display a progress bar during the spatial loop. Default is True.

    notebook : bool, optional
        If True, use the Jupyter notebook version of the progress bar.
    

    Returns
    -------
    None
        The resulting datacube is written directly to disk.

    Outputs
    -------
    Primary HDU
        Cropped or spatially resampled flux datacube.

    ``Error_cube``
        Corresponding uncertainty datacube.

    ``BADPIXELMASK``
        Integer mask containing 1 for pixels mapped inside ``file1`` and
        0 for pixels outside its spatial footprint.

    Notes
    -----
    When ``spsample_copy=False``, the FoV defined by ``file0`` is
    transformed into the pixel coordinate system of ``file1`` and the
    corresponding rectangular section of ``file1`` is extracted. The
    output therefore retains the spatial sampling and orientation of
    ``file1``.

    When ``spsample_copy=True``, an output spatial grid is constructed
    from the WCS of ``file0``. Each output pixel is transformed to sky
    coordinates and subsequently into the pixel system of ``file1``.
    Flux and uncertainty values are evaluated using
    ``tools.map_interpolB``.

    ``fac_sizeX`` and ``fac_sizeY`` modify the dimensions of the
    reference grid while keeping its central position fixed.

    ``dx`` and ``dy`` shift the astrometric centre of the reference FoV
    in units of reference-cube pixels before either cropping or
    resampling is performed.

    Examples
    --------
    Crop ``file1`` using the FoV of ``file0`` while preserving the
    spatial sampling of ``file1``::

        >>> crop_cube(
        ...     'reference.fits.gz',
        ...     'input.fits.gz',
        ...     'crop.fits.gz'
        ... )

    Resample the output onto the spatial grid of the reference cube::

        >>> crop_cube(
        ...     'reference.fits.gz',
        ...     'input.fits.gz',
        ...     'crop.fits.gz',
        ...     spsample_copy=True
        ... )

    Use a FoV 1.5 times larger than the reference FoV::

        >>> crop_cube(
        ...     'reference.fits.gz',
        ...     'input.fits.gz',
        ...     'crop.fits.gz',
        ...     fac_sizeX=1.5,
        ...     fac_sizeY=1.5
        ... )

    Apply an astrometric displacement of two pixels in x and -1 pixel
    in y::

        >>> crop_cube(
        ...     'reference.fits.gz',
        ...     'input.fits.gz',
        ...     'crop.fits.gz',
        ...     dx=2,
        ...     dy=-1
        ... )

    See Also
    --------
    crop_image
        Reproject external images onto the spatial grid of a datacube.

    coad_cube
        Co-add datacubes from multiple spectral bands.
    """

    # ---------------------------------------------------------
    # Check input parameters
    # ---------------------------------------------------------

    if fac_sizeX <= 0 or fac_sizeY <= 0:
        raise ValueError(
            "fac_sizeX and fac_sizeY must be greater than zero."
        )

    # ---------------------------------------------------------
    # Read reference cube
    # ---------------------------------------------------------
    
    print("Reading reference cube:", file0)
    try:
        cube0, hdr0 = fits.getdata(
            file0, 0, header=True)
    except:
        cube0, hdr0 = fits.getdata(
            file0, 1, header=True)

    if cube0.ndim != 3:
        raise ValueError(
            "The primary HDU of file0 must contain a 3D datacube."
        )

    nz0, nx0, ny0 = cube0.shape

    wcs0 = WCS(hdr0).celestial

    # Determine input-reference pixel area.
    dxt = np.sqrt(hdr0['CD1_1']**2.0 +
                 hdr0['CD1_2']**2.0) * 3600.0
    dyt = np.sqrt(hdr0['CD2_1']**2.0 +
                 hdr0['CD2_2']**2.0) * 3600.0

    A0 = dxt * dyt

    # ---------------------------------------------------------
    # Read input cube
    # ---------------------------------------------------------

    print("Reading input cube:", file1)
    
    try:
        cube1, hdr1 = fits.getdata(
            file1, 0, header=True
        )

        cube1E = fits.getdata(
            file1, 1, header=False
        )
    except:
        cube1, hdr1 = fits.getdata(
            file1, 1, header=True
        )

        cube1E = fits.getdata(
            file1, 2, header=False
        )

    if cube1.ndim != 3:
        raise ValueError(
            "The primary HDU of file1 must contain a 3D datacube."
        )

    if cube1E.shape != cube1.shape:
        raise ValueError(
            "The uncertainty cube must have the same shape "
            "as the flux cube."
        )

    nz1, nx1, ny1 = cube1.shape

    wcs1 = WCS(hdr1).celestial

    # Determine the reference pixel scale and pixel area.
    dxt = np.sqrt(hdr1['CD1_1']**2.0 +
                 hdr1['CD1_2']**2.0) * 3600.0
    dyt = np.sqrt(hdr1['CD2_1']**2.0 +
                 hdr1['CD2_2']**2.0) * 3600.0
    A1 = dxt * dyt

    cube1 = cube1 / A1 * A0
    cube1E = cube1E / A1 * A0
    

    # ---------------------------------------------------------
    # Define reference FoV
    # ---------------------------------------------------------
    #
    # FITS x corresponds to the third NumPy dimension (ny).
    # FITS y corresponds to the second NumPy dimension (nx).
    #
    # The geometrical centre is shifted by dx/dy reference pixels.
    # ---------------------------------------------------------

    xcen = (ny0 - 1) / 2.0 + dx
    ycen = (nx0 - 1) / 2.0 + dy

    size_x = ny0 * fac_sizeX
    size_y = nx0 * fac_sizeY

    # Pixel edges of the desired FoV.
    x_min_ref = xcen - size_x / 2.0
    x_max_ref = xcen + size_x / 2.0

    y_min_ref = ycen - size_y / 2.0
    y_max_ref = ycen + size_y / 2.0

    # Four corners are used to support rotated WCS solutions.
    corners_x = np.array([
        x_min_ref,
        x_max_ref,
        x_max_ref,
        x_min_ref
    ])

    corners_y = np.array([
        y_min_ref,
        y_min_ref,
        y_max_ref,
        y_max_ref
    ])

    sky_corners = pixel_to_skycoord(
        corners_x,
        corners_y,
        wcs0
    )

    # =========================================================
    # MODE 1: preserve spatial sampling of file1
    # =========================================================

    if not spsample_copy:

        xpos, ypos = skycoord_to_pixel(
            sky_corners,
            wcs1
        )

        valid = (
            np.isfinite(xpos)
            & np.isfinite(ypos)
        )

        if not np.any(valid):
            raise ValueError(
                "The requested FoV cannot be transformed "
                "into the WCS of file1."
            )

        xpos = xpos[valid]
        ypos = ypos[valid]

        # Spatial limits in file1.
        ymin = int(np.floor(np.min(xpos)))
        ymax = int(np.ceil(np.max(xpos))) + 1

        xmin = int(np.floor(np.min(ypos)))
        xmax = int(np.ceil(np.max(ypos))) + 1

        # Check whether the two FoVs overlap.
        if (
            xmax <= 0
            or ymax <= 0
            or xmin >= nx1
            or ymin >= ny1
        ):
            raise ValueError(
                "The requested FoV does not overlap file1."
            )

        # Restrict to file1 boundaries.
        xmin = max(0, xmin)
        xmax = min(nx1, xmax)

        ymin = max(0, ymin)
        ymax = min(ny1, ymax)

        print(
            "Cropping limits in file1: "
            "x=[{}, {}], y=[{}, {}]".format(
                xmin,
                xmax,
                ymin,
                ymax
            )
        )

        # Direct NumPy extraction.
        cube_out = cube1[
            :,
            xmin:xmax,
            ymin:ymax
        ].copy()

        cube_outE = cube1E[
            :,
            xmin:xmax,
            ymin:ymax
        ].copy()

        cube_outB = np.ones(
            cube_out.shape,
            dtype=int
        )

        # Output header follows file1.
        hdr_out = hdr1.copy()

        if 'CRPIX1' in hdr_out:
            hdr_out['CRPIX1'] -= ymin

        if 'CRPIX2' in hdr_out:
            hdr_out['CRPIX2'] -= xmin

    # =========================================================
    # MODE 2: copy spatial sampling of file0
    # =========================================================

    else:

        # Number of output pixels.
        ny_out = max(
            1,
            int(np.round(ny0 * fac_sizeX))
        )

        nx_out = max(
            1,
            int(np.round(nx0 * fac_sizeY))
        )

        print(
            "Resampling cube to reference spatial grid: "
            "{} x {} pixels".format(
                nx_out,
                ny_out
            )
        )

        # -----------------------------------------------------
        # Construct output header
        # -----------------------------------------------------
        #
        # Spectral WCS comes from file1.
        # Spatial WCS comes from file0.
        # -----------------------------------------------------

        hdr_out = hdr1.copy()

        # Copy the celestial WCS keywords from file0.
        spatial_keys = [
            'CTYPE1', 'CTYPE2',
            'CUNIT1', 'CUNIT2',
            'CRVAL1', 'CRVAL2',
            'CRPIX1', 'CRPIX2',
            'CDELT1', 'CDELT2',
            'CD1_1', 'CD1_2',
            'CD2_1', 'CD2_2',
            'PC1_1', 'PC1_2',
            'PC2_1', 'PC2_2'
        ]

        for key in spatial_keys:

            if key in hdr0:
                hdr_out[key] = hdr0[key]

        # -----------------------------------------------------
        # Adjust CRPIX for enlarged FoV and astrometric shift
        # -----------------------------------------------------
        #
        # The new array is centred on the shifted reference
        # position.
        # -----------------------------------------------------

        xcen_out = (ny_out - 1) / 2.0
        ycen_out = (nx_out - 1) / 2.0

        # Sky coordinate of the shifted reference centre.
        sky_centre = pixel_to_skycoord(
            xcen - dx,
            ycen - dy,
            wcs0
        )

        # Build an intermediate celestial WCS with the same
        # scale/orientation as file0.
        wcs_out = WCS(hdr_out).celestial

        # Determine where the desired centre currently lies.
        x_tmp, y_tmp = skycoord_to_pixel(
            sky_centre,
            wcs_out
        )

        # Shift CRPIX so that sky_centre falls exactly at the
        # geometrical centre of the output cube.
        hdr_out['CRPIX1'] += (
            xcen_out - x_tmp
        )

        hdr_out['CRPIX2'] += (
            ycen_out - y_tmp
        )

        # Reconstruct WCS after changing CRPIX.
        wcs_out = WCS(hdr_out).celestial

        # -----------------------------------------------------
        # Allocate output arrays
        # -----------------------------------------------------

        cube_out = np.full(
            (nz1, nx_out, ny_out),
            np.nan,
            dtype=float
        )

        cube_outE = np.full(
            (nz1, nx_out, ny_out),
            np.nan,
            dtype=float
        )

        cube_outB = np.zeros(
            (nz1, nx_out, ny_out),
            dtype=int
        )

        # -----------------------------------------------------
        # Spatial resampling
        # -----------------------------------------------------

        if pbars:
            if notebook:
                pbar=tqdm(total=nx_out)
            else:     
                pbar=tqdmT(total=nx_out) 

        for i in range(nx_out):

            for j in range(ny_out):

                # Output pixel -> sky.
                sky = pixel_to_skycoord(
                    j + dx,
                    i + dy,
                    wcs_out
                )

                # Sky -> input cube pixel.
                xpos, ypos = skycoord_to_pixel(
                    sky,
                    wcs1
                )

                # Require enough space for interpolation.
                if (
                    xpos >= 0
                    and xpos < ny1
                    and ypos >= 0
                    and ypos < nx1
                ):

                    # Interpolate every wavelength plane.

                    cube_out[:, i, j] = \
                        tools.cube_interpolB(
                            cube1,
                            ypos,
                            xpos
                        )

                    cube_outE[:, i, j] = \
                        tools.cube_interpolB(
                            cube1E,
                            ypos,
                            xpos
                        )

                    cube_outB[:, i, j] = 1
            if pbars:
                pbar.update(1)
        if pbars:
            pbar.close()          

    # ---------------------------------------------------------
    # Construct FITS HDUs
    # ---------------------------------------------------------

    h1 = fits.PrimaryHDU(
        cube_out,
        header=hdr_out.copy()
    )

    h2 = fits.ImageHDU(
        cube_outE,
        header=hdr_out.copy(),
        name='Error_cube'
    )

    h3 = fits.ImageHDU(
        cube_outB,
        header=hdr_out.copy(),
        name='BADPIXELMASK'
    )

    h2.header['EXTNAME'] = 'Error_cube'
    h3.header['EXTNAME'] = 'BADPIXELMASK'

    # ---------------------------------------------------------
    # Output filename
    # ---------------------------------------------------------

    output_file = str(file2)

    if not (
        output_file.endswith('.fits')
        or output_file.endswith('.fits.gz')
    ):
        output_file += '.fits'

    # ---------------------------------------------------------
    # Write output
    # ---------------------------------------------------------

    hlist = fits.HDUList([
        h1,
        h2,
        h3
    ])

    hlist.writeto(
        output_file,
        overwrite=True
    )

    tools.sycall(
        'gzip -f ' + output_file
    )

    print("Output cube:", output_file)


def extract_cube_region(file0, file2, file1=None, reg_file='Reg.reg',
                        dir_reg='./', mask=False, pbar=True,
                        notebook=True):
    """
    Extract a spatial subcube from an IFU datacube using a DS9 region.

    The function reads a rectangular DS9 region, converts its celestial
    coordinates to the pixel coordinate system of the input datacube,
    and extracts the corresponding spatial region while preserving the
    complete spectral axis.

    The flux, uncertainty, and bad-pixel-mask cubes are written to a new
    FITS file. The spatial WCS reference pixels are updated so that the
    celestial coordinates of the extracted cube remain consistent with
    those of the original datacube.

    Optionally, a second reference cube can be used as a spatial mask.
    Pixels in the extracted region that do not overlap valid data in the
    reference cube are set to NaN and flagged in the bad-pixel mask.

    Parameters
    ----------
    file0 : str or path-like
        Input FITS datacube. The primary HDU must contain the flux cube
        with dimensions ``(wavelength, x, y)``.

        Extension 1 is expected to contain the corresponding uncertainty
        cube with the same dimensions as the primary HDU.

    file2 : str or path-like
        Output filename. The ``.fits`` or ``.fits.gz`` extension may be
        included or omitted. The output file is written directly as a
        FITS file.

    file1 : str or path-like, optional
        Reference datacube used to define a spatial validity mask when
        ``mask=True``. Its primary HDU is collapsed along the spectral
        axis and reprojected through the WCS to determine whether each
        output spatial pixel overlaps valid reference data.

        If ``mask=False``, this parameter is ignored. Default is None.

    reg_file : str, optional
        Name of the DS9 region file defining the spatial region to
        extract. The first region in the file is used. The current
        implementation expects a rectangular/box aperture for which
        ``get_apertures`` provides the central coordinates and spatial
        dimensions. Default is ``'Reg.reg'``.

    dir_reg : str or path-like, optional
        Directory containing ``reg_file``. Default is ``'./'``.

    mask : bool, optional
        If True, use ``file1`` to determine which spatial pixels overlap
        valid data in the reference cube. Pixels outside the reference
        data footprint are set to NaN and flagged in the bad-pixel mask.
        Default is False.

    pbar : bool, optional
        If True, display a progress bar during the spatial extraction.
        Default is True.

    notebook : bool, optional
        If True, use the Jupyter notebook version of the tqdm progress
        bar. If False, use the terminal version. Default is True.

    Returns
    -------
    None
        The function writes the extracted datacube directly to disk.

    Outputs
    -------
    Primary HDU
        Extracted flux datacube.

    ``Error_cube``
        Extracted uncertainty datacube.

    ``BADPIXELMASK``
        Integer mask with a value of 1 for retained pixels and 0 for
        pixels rejected by the optional reference-cube mask.

    Notes
    -----
    The celestial limits of the extraction are calculated from the
    centre and dimensions of the first aperture returned by
    :func:`CubeGen.tools.tools.get_apertures`.

    The region limits are transformed from celestial coordinates to
    pixels using the celestial WCS of the input cube.

    After extraction, the FITS reference pixels are shifted according
    to the origin of the extracted region,

    .. math::

        \\mathrm{CRPIX1}_{new}
        =
        \\mathrm{CRPIX1}_{old} - y_{min}

    and

    .. math::

        \\mathrm{CRPIX2}_{new}
        =
        \\mathrm{CRPIX2}_{old} - x_{min}.

    This preserves the original celestial coordinate system in the
    extracted datacube.

    When ``mask=True``, every spatial pixel is transformed from the
    input-cube WCS to the reference-cube WCS. The collapsed reference
    image is evaluated using
    :func:`CubeGen.tools.tools.map_interpolB`.

    Examples
    --------
    Extract a region without applying an external spatial mask::

        >>> extract_cube_region(
        ...     'galaxy_cube.fits.gz',
        ...     'galaxy_region.fits.gz',
        ...     reg_file='Reg.reg'
        ... )

    Extract a region and restrict it to the footprint of another cube::

        >>> extract_cube_region(
        ...     'galaxy_cube.fits.gz',
        ...     'galaxy_region.fits.gz',
        ...     file1='reference_cube.fits.gz',
        ...     reg_file='Reg.reg',
        ...     mask=True
        ... )

    See Also
    --------
    crop_image
        Reproject external imaging onto the spatial grid of an IFU cube.

    coad_cube
        Co-add reconstructed datacubes from multiple spectral bands.
    """

    # ---------------------------------------------------------
    # Read DS9 aperture
    # ---------------------------------------------------------

    reg_path = dir_reg + reg_file

    ra_R, dec_R, rad_R, l1_R, l2_R, th_R, color, names, typ = \
        tools.get_apertures(reg_path)

    if len(ra_R) == 0:
        raise ValueError(
            "No valid aperture was found in {}".format(reg_path)
        )

    # Use the first aperture.
    ra = ra_R[0]
    dec = dec_R[0]
    l1 = l1_R[0]
    l2 = l2_R[0]

    # ---------------------------------------------------------
    # Read input cube
    # ---------------------------------------------------------

    print("Reading cube:", file0)

    pdl_cube0, hdr0 = fits.getdata(
        file0, 0, header=True
    )

    pdl_cube0E = fits.getdata(
        file0, 1, header=False
    )

    if pdl_cube0.ndim != 3:
        raise ValueError(
            "The primary HDU of file0 must contain a 3D datacube."
        )

    if pdl_cube0E.shape != pdl_cube0.shape:
        raise ValueError(
            "The uncertainty cube must have the same shape "
            "as the flux cube."
        )

    nz0, nx0, ny0 = pdl_cube0.shape

    # ---------------------------------------------------------
    # Convert region centre and limits to celestial coordinates
    # ---------------------------------------------------------

    # Local imports required by this function.
    from astropy.coordinates import SkyCoord, FK5
    from astropy import units as u

    sky_centre = SkyCoord(
        ra + ' ' + dec,
        frame=FK5,
        unit=(u.hourangle, u.deg)
    )

    ra_deg = sky_centre.ra.deg
    dec_deg = sky_centre.dec.deg

    # Region dimensions are given in arcseconds.
    ra1 = ra_deg - l1 / 2.0 / 3600.0
    ra2 = ra_deg + l1 / 2.0 / 3600.0

    dec1 = dec_deg - l2 / 2.0 / 3600.0
    dec2 = dec_deg + l2 / 2.0 / 3600.0

    sky00 = SkyCoord(
        ra1, dec1,
        frame=FK5,
        unit=(u.deg, u.deg)
    )

    sky11 = SkyCoord(
        ra2, dec2,
        frame=FK5,
        unit=(u.deg, u.deg)
    )

    # ---------------------------------------------------------
    # Transform region limits to input-cube pixels
    # ---------------------------------------------------------

    wcs0 = WCS(hdr0).celestial

    ypos00, xpos00 = skycoord_to_pixel(
        sky00, wcs0
    )

    ypos11, xpos11 = skycoord_to_pixel(
        sky11, wcs0
    )

    xpos00 = int(np.round(xpos00))
    ypos00 = int(np.round(ypos00))

    xpos11 = int(np.round(xpos11))
    ypos11 = int(np.round(ypos11))

    # Ensure increasing array limits.
    xmin = min(xpos00, xpos11)
    xmax = max(xpos00, xpos11)

    ymin = min(ypos00, ypos11)
    ymax = max(ypos00, ypos11)

    # Restrict extraction to the input-cube boundaries.
    xmin = max(0, xmin)
    xmax = min(nx0, xmax)

    ymin = max(0, ymin)
    ymax = min(ny0, ymax)

    if xmax <= xmin or ymax <= ymin:
        raise ValueError(
            "The requested region does not overlap the input cube."
        )

    print(
        "Extraction limits: "
        "x=[{}, {}], y=[{}, {}]".format(
            xmin, xmax, ymin, ymax
        )
    )

    # ---------------------------------------------------------
    # Allocate output cubes
    # ---------------------------------------------------------

    nx2 = xmax - xmin
    ny2 = ymax - ymin

    seg = np.full(
        (nz0, nx2, ny2),
        np.nan,
        dtype=float
    )

    seg_e = np.full(
        (nz0, nx2, ny2),
        np.nan,
        dtype=float
    )

    seg_B = np.ones(
        (nz0, nx2, ny2),
        dtype=int
    )

    # ---------------------------------------------------------
    # Optional reference-cube mask
    # ---------------------------------------------------------

    if mask:

        if file1 is None:
            raise ValueError(
                "file1 must be provided when mask=True."
            )

        print("Reading reference cube:", file1)

        pdl_cube1, hdr1 = fits.getdata(
            file1, 0, header=True
        )

        if pdl_cube1.ndim == 3:
            map1 = np.nansum(
                pdl_cube1,
                axis=0
            )
        elif pdl_cube1.ndim == 2:
            map1 = np.copy(pdl_cube1)
        else:
            raise ValueError(
                "file1 must contain either a 2D image "
                "or a 3D datacube."
            )

        nx1, ny1 = map1.shape

        wcs1 = WCS(hdr1).celestial

    # ---------------------------------------------------------
    # Extract region
    # ---------------------------------------------------------

    iterator = range(xmin, xmax)

    if pbar:

        if notebook:
            iterator = tqdm(
                iterator,
                total=nx2,
                desc='Extracting cube'
            )
        else:
            iterator = tqdmT(
                iterator,
                total=nx2,
                desc='Extracting cube'
            )

    for i in iterator:

        for j in range(ymin, ymax):

            valid = True

            # -------------------------------------------------
            # Check reference-cube footprint
            # -------------------------------------------------

            if mask:

                sky = pixel_to_skycoord(
                    j, i, wcs0
                )

                xpos_ref, ypos_ref = skycoord_to_pixel(
                    sky, wcs1
                )

                if (
                    xpos_ref < 0
                    or xpos_ref >= ny1
                    or ypos_ref < 0
                    or ypos_ref >= nx1
                ):

                    valid = False

                else:

                    val = tools.map_interpolB(
                        map1,
                        ypos_ref,
                        xpos_ref
                    )

                    if (
                        not np.isfinite(val)
                        or val == 0
                    ):
                        valid = False

            # -------------------------------------------------
            # Copy flux and uncertainty
            # -------------------------------------------------

            ii = i - xmin
            jj = j - ymin

            if valid:

                seg[:, ii, jj] = \
                    pdl_cube0[:, i, j]

                seg_e[:, ii, jj] = \
                    pdl_cube0E[:, i, j]

            else:

                seg[:, ii, jj] = np.nan
                seg_e[:, ii, jj] = np.nan
                seg_B[:, ii, jj] = 0

    # ---------------------------------------------------------
    # Construct FITS output
    # ---------------------------------------------------------

    h1 = fits.PrimaryHDU(seg)
    h2 = fits.ImageHDU(seg_e)
    h3 = fits.ImageHDU(seg_B)

    # Copy original header.
    keys = list(hdr0.keys())

    for hdu in (h1, h2, h3):

        hdr = hdu.header

        for key in keys:

            try:
                hdr[key] = hdr0[key]
                hdr.comments[key] = hdr0.comments[key]
            except:
                pass

        # Shift reference pixels to the new spatial origin.
        hdr['CRPIX1'] = hdr0['CRPIX1'] - ymin
        hdr['CRPIX2'] = hdr0['CRPIX2'] - xmin

    h2.header['EXTNAME'] = 'Error_cube'
    h3.header['EXTNAME'] = 'BADPIXELMASK'

    # ---------------------------------------------------------
    # Output filename
    # ---------------------------------------------------------

    output_file = str(file2)

    if output_file.endswith('.fits.gz'):
        pass
    elif output_file.endswith('.fits'):
        pass
    else:
        output_file += '.fits'

    # ---------------------------------------------------------
    # Write cube
    # ---------------------------------------------------------

    hlist = fits.HDUList(
        [h1, h2, h3]
    )

    hlist.writeto(
        output_file,
        overwrite=True
    )

    print("Output cube:", output_file)