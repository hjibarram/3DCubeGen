from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import CubeGen.megaratools.megtools as mtools
import CubeGen.tools.tools as tools
import numpy as np


def astromatch(file0,file1,sig=2):
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
    nz0,nx0,ny0=spec0.shape
    
    [spec1, hdr1]=fits.getdata(file1, 0, header=True)
    nz1,nx1,ny1=spec1.shape

    # Collapse the cubes along the spectral axis.
    try:
        map0=np.nansum(spec0,axis=0)
    except:
        map0=np.copy(spec0)
    try:
        map1=np.nansum(spec1,axis=0)
    except:
        map1=np.copy(spec1)
    print(file0)
    print(file1)

    # Determine the PSF centroid in each reconstructed image.
    x0,y0,ds_m0,psf0,model0=mtools.evaluate_2dPSF(map0,model=True,sig=sig)
    x1,y1,ds_m1,psf1,model1=mtools.evaluate_2dPSF(map1,model=True,sig=sig)

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


def crop_image(names, cube, dir1='.', dir2='.', apt='_gri'):
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

                sky1 = pixel_to_skycoord(i, j, wcs1)

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
        cube.replace('.fits.gz', apt) + '.fits',
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

    from PIL import Image

    im = Image.fromarray(rgb_cube)

    im.save(
        cube.replace('.fits.gz', apt) + '.jpeg',
        quality=100
    )