import numpy as np
from hoki import load
import glob
import ppxf.ppxf_util as util
from scipy import ndimage

def get_BPASS_library(filepath, gal_lamrange):
    """
    Loads the first BPASS model into a numpy array, gets the library filenames list,
    the wavelength range of the SSPs and the list of ages

    Parameters
    -----------
    filepath : str
        The filenames of the models to be loaded

    gal_lamrange : :obj:'~numpy.ndarray'
        The first and last wavelengths of the galaxy data going to be fit

    Returns
    -------
    ssp_lamrange : :obj:'~numpy.ndarray'
        Wavelength range of the SSPs

    templates : :obj:'~numpy.ndarray'
        Array of BPASS templates ready to be input into ppxf, needs to be
        TEMPLATES[nPixels, nAge, nMetal]

    ssp_ages : list
        List of ages of the templates

    ssp_metals : list
        List of metals of the templates
    """
    #first we need to make a list of the filenames
    ssp_lib = glob.glob(filepath)

    #open the first file to get the wavelength vector and shape for the rest of the array
    #this creates a pandas dataframe
    ssp_data = load.model_output(ssp_lib[0])

    #get the log ages
    ssp_ages = list(ssp_data.columns)[1:]

    #make into a numpy array
    ssp_data = ssp_data.to_numpy()

    ssp_lamdas = ssp_data[:,0]
    ssp_data = ssp_data[:,1:]

    #mask the wavelengths 3000A-6000A, since the wavelength range of KCWI is 350-560nm
    lamda_mask = (ssp_lamdas>gal_lamrange[0]-200)&(ssp_lamdas<gal_lamrange[1]+200)
    ssp_lamdas = ssp_lamdas[lamda_mask]
    ssp_data = ssp_data[lamda_mask,:]

    #get the lam_range
    ssp_lamrange = np.array([ssp_lamdas[0], ssp_lamdas[-1]])

    #create empty templates array to add templates to in shape [nPix, nAge, nMetal]
    templates = np.empty([ssp_data.shape[0], ssp_data.shape[1], len(ssp_lib)])

    #create empty list to add the metals to
    ssp_metals = []

    for j, file in enumerate(ssp_lib):
        #load the file
        ssp_data = load.model_output(file)
        #convert to numpy array and ignore the wavelength vector
        ssp_data = ssp_data.to_numpy()[:,1:]
        #add to templates array
        templates[:,:,j] = ssp_data[lamda_mask,:]
        #add metals to metals list
        ssp_metals.append(file[-11:-7])

    print('Templates shape: ', templates.shape)

    return ssp_lamrange, templates, ssp_ages, ssp_metals




def smoothing(fwhm1, fwhm2, input_spectrum, cdelt):
    """
    Smooth either the template to the spectrum, or the spectrum to the template
    FWHM. fwhm1 must always be larger than fwhm2.  Usually want to be smoothing
    template to the spectrum - should have higher resolution templates than data.

    Parameters
    ----------
    fwhm1 : float or :obj: '~numpy.ndarray'
        the first and larger FWHM value, should be the one that you want to
        smooth to (the data)

    fwhm2 : float or :obj: '~numpy.ndarray'
        the second and smaller FWHM value, the one you are smoothing (from the
        templates)

    input_spectrum : :obj: '~numpy.ndarray'
        the spectrum to smooth (usually the template)

    cdelt : float
        the wavelength Angstroms per pixel from the templates (might be CDELT1
        or CD3_3 in a fits header)

    Returns
    -------
    input_spectrum : :obj: '~numpy.ndarray'
        the spectrum smoothed to fwhm1
    """
    #find the difference in the fwhm
    fwhm_diff = np.sqrt(fwhm1**2 - fwhm2**2)
    #calculate the sigma difference in pixels
    sigma = fwhm_diff/2.355/cdelt

    #use the difference to smooth the input
    #if the fwhm1 is the same for all pixels, use the scipy convolution,
    #otherwise use the one from ppxf, which can handle convolving a spectrum
    #by a Gaussian with a different sigma for every pixel
    if np.isscalar(fwhm1):
        input_spectrum = ndimage.gaussian_filter1d(input_spectrum, sigma)
    else:
        input_spectrum = gaussian_filter1d(input_spectrum, sigma)

    return input_spectrum

def gaussian_filter1d(spec, sig):
    """
    Convolve a spectrum by a Gaussian with different sigma for every pixel.
    If all sigma are the same this routine produces the same output as
    scipy.ndimage.gaussian_filter1d, except for the border treatment.
    Here the first/last p pixels are filled with zeros.
    When creating a template library for SDSS data, this implementation
    is 60x faster than a naive for loop over pixels.

    :param spec: vector with the spectrum to convolve
    :param sig: vector of sigma values (in pixels) for every pixel
    :return: spec convolved with a Gaussian with dispersion sig

    """
    sig = sig.clip(0.01)  # forces zero sigmas to have 0.01 pixels
    p = int(np.ceil(np.max(3*sig)))
    m = 2*p + 1  # kernel size
    x2 = np.linspace(-p, p, m)**2

    n = spec.size
    a = np.zeros((m, n))
    for j in range(m):   # Loop over the small size of the kernel
        a[j, p:-p] = spec[j:n-m+j+1]

    gau = np.exp(-x2[:, None]/(2*sig**2))
    gau /= np.sum(gau, 0)[None, :]  # Normalize kernel

    conv_spectrum = np.sum(a*gau, 0)

    return conv_spectrum

def prep_BPASS_models(ssp_templates, ssp_lamrange, gal_velscale, lamrange_gal, z, fwhm_gal=1.7, fwhm_temp=2.0, cdelt_temp=1.0, velscale_ratio=1, em_lines=False, fwhm_emlines=2.0, vacuum=True, extra_em_lines=False, tie_balmer=True):
    """
    Smooths, Log rebins and normalises the templates

    Parameters
    ----------
    ssp_templates : :obj:'~numpy.ndarray'
        The templates in a numpy array

    ssp_lamrange : :obj:'~numpy.ndarray'
        Wavelength range of the SSPs

    gal_velscale : float
        the velocity scale of the data

    lamrange_gal : :obj:'~numpy.ndarray'
        the wavelength range of the data

    z : float
        the redshift of the galaxy

    fwhm_gal : float
        the fwhm of the data

    fwhm_temp : float
        the fwhm of the templates

    cdelt_temp : float
        the delta wavelength of the templates in Angstroms/pixel

    velscale_ratio : float
        the number by which to divide the galaxy velscale when log rebinning the
        SSPs (default=1, gives the same spectral sampling for the templates and
        the galaxy data. 2 would give templates twice sampled compared to galaxy
        data.)

    em_lines : boolean
        whether to include emission lines in the templates

    fwhm_emlines : float
        the fwhm of the emission lines added to the templates

    vacuum : boolean
        whether to make the wavelengths in vacuum or air wavelengths
        (default=True, default in ppxf is False.)

    extra_em_lines : bool
        set to True to include extra emission lines often found in
        KCWI data (OII 4317, [OIII]4363, OII4414 and [NeIII]3868).
        (default=False)

    tie_balmer : bool
        ties the Balmer lines according to a theoretical decrement
        (case B recombination T=1e4 K, n=100 cm^-3) (default=True)

    Returns
    -------
    templates : :obj: '~numpy.ndarray'

    ssp_logLam : :obj:'~numpy.ndarray'
        the logarithmically rebinned wavelength vector for the templates and
        emission lines (if included)

    Notes
    -----
    If emission lines are included in the fit also Returns:
        gas_names : the names of the lines included
        line_wave : the wavelength of the lines included in vacuum wavelengths
        component : assigns which templates belong to which components (stellar,
        emission lines, Balmer lines)
        gas_component : vector, True for gas templates
        temp_dim : the original dimensions of the templates without the emission
        lines
    """
    #apply log_rebin to the first template, and sample at a velocity scale which defaults to the same as that for the spectra
    log_ssp, ssp_loglam, ssp_velscale = util.log_rebin(ssp_lamrange, ssp_templates[:,0,0], velscale=gal_velscale/velscale_ratio)
    print('rebin', log_ssp)
    #create an array for the templates
    templates = np.empty([log_ssp.size, ssp_templates.shape[1], ssp_templates.shape[2]])
    print('temp', templates)
    #iterate through all of the templates and add to the empty array
    for i in np.arange(templates.shape[1]):
        for j in np.arange(templates.shape[2]):
            #smooth the template if the template fwhm is less than the galaxy fwhm
            if fwhm_gal > fwhm_temp:
                ssp = smoothing(fwhm_gal, fwhm_temp, ssp_templates[:,j], cdelt_temp)
            else:
                ssp = ssp_templates[:,j]
            #log-rebin the template
            ssp = util.log_rebin(ssp_lamrange, ssp, velscale=gal_velscale/velscale_ratio)[0]
            print('ssp', ssp)
            print('temp', templates)
            #take the median and add to templates
            templates[:,i,j] = ssp

    #divide all the templates by the median
    templates = templates/np.nanmedian(templates)
    templates = templates.reshape([templates.shape[0], -1])
    return templates, ssp_loglam

