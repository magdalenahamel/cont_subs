'''
script to run only one spaxel on ppxf'''


import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import ppxf.ppxf_util as util
from load_model import get_BPASS_library
from load_model import prep_BPASS_models
from load_model import smoothing
import glob
from hoki import load



#Load Data
pix = [35, 62]
# Why doesnt work open the file from Data directory
#data_folder = '/Users/magdalenahamel/Desktop/PhD/Data/ncg1569/'
#file = '/Users/magdalenahamel/Desktop/PhD/Data/ncg1569/ngc1569_pointing_2_red_metacube_mw.fits'

fwhm_gal=1.7
fwhm_temp=2.0
file = 'ngc1569_pointing_2_red_metacube_mw.fits'
hdu = fits.open(file)
head = hdu[0].header
cube = hdu[0].data

spec = cube[:, pix[0], pix[1]]
npix = cube.shape[0]
wave = head['CRVAL3'] + head['CDELT3']*np.arange(npix)
noise = np.full_like(spec, 0.0047)


#SN array
cont_mask = (wave>4600) & (wave<4800)
sn_array = abs(spec[cont_mask]/noise[cont_mask])
#for cube use: sn_array = abs(np.nanmedian(data_flat[cont_mask,:]/noise_flat[cont_mask,:], axis=0))


#if the data FWHM is less than the template FWHM, we can choose to smooth the
#galaxy data (this is not a good thing though)
if fwhm_gal < fwhm_temp:
    print('WARNING:     SMOOTHING DATA TO LOWER RESOLUTION')
    spec = smoothing(fwhm_temp, fwhm_gal, spec, cdelt=0.5)
    noise = smoothing(fwhm_temp, fwhm_gal, noise, cdelt=0.5)

#Mask data so it doesnt cover wavelenghts thata are not in models
#not necessary for this case

#Logbin and normalize spectra
spec_logspec = np.empty_like(spec)
noise_logspace = np.empty_like(noise)
gal_velscale = np.empty_like(spec)

wave_range = np.array([wave[0], wave[-1]])

spec_log, wave_log, velscale = util.log_rebin(wave, spec)
noise_log, noise_wave_log, noise_velscale = util.log_rebin(wave, noise)

#normalise the spectra and noise
spec_norm = np.nanmedian(spec_log, axis=0)
spec_log = spec_log/spec_norm
noise_norm = np.nanmedian(noise_log, axis=0)
noise_log = noise_log/noise_norm

###Load the models
filepath = '/Users/magdalenahamel/Desktop/PhD/Data/models/*'
ssp_lib = glob.glob(filepath) #list with the filepaths
#load one to get wavelenght
ssp_data = load.model_output(ssp_lib[0])
# get the log ages
ssp_ages = list(ssp_data.columns)[1:]

#make into a numpy array
ssp_data = ssp_data.to_numpy()
ssp_lamdas = ssp_data[:,0]
ssp_data = ssp_data[:,1:]

# mask the wavelengths 3000A-6000A, since the wavelength range of KCWI is 350-560nm
lamda_mask = (ssp_lamdas > wave[0] - 200) & (ssp_lamdas < wave[-1] + 200)
ssp_lamdas = ssp_lamdas[lamda_mask]
ssp_data = ssp_data[lamda_mask, :]

# get the lam_range
ssp_lamrange = np.array([ssp_lamdas[0], ssp_lamdas[-1]])

# create empty templates array to add templates to in shape [nPix, nAge, nMetal]
templates = np.empty([ssp_data.shape[0], ssp_data.shape[1], len(ssp_lib)])

# create empty list to add the metals to
ssp_metals = []

for j, file in enumerate(ssp_lib):
    # load the file
    ssp_data = load.model_output(file)
    # convert to numpy array and ignore the wavelength vector
    ssp_data = ssp_data.to_numpy()[:, 1:]
    # add to templates array
    templates[:, :, j] = ssp_data[lamda_mask, :]
    # add metals to metals list
    ssp_metals.append(file[-11:-7])

print('Templates shape: ', templates.shape)
ssp_lamrange, templates, ssp_ages, ssp_metals

## prepare templates

templates, ssp_logLam, gas_names, line_wave, component, gas_component, temp_dim

