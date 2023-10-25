import matplotlib.pyplot as plt
import numpy as np
from astropy.table import QTable
from photutils.datasets import (make_gaussian_sources_image,
                                make_noise_image)
import matplotlib
matplotlib.use('TkAgg')
from astropy.io import fits
from libraries.utils import Utils

dir = "D:\\July2022\\2022-08-21\\test_star\\raw\\transmission\\"
shape = (2048, 2048)

# set up the poisiton and flux parameters to tweak
x_var = np.random.normal(900, 2, 3600)
y_var = np.random.normal(900, 2, 3600)

x_star = np.random.normal(1300, 2, 3600)
y_star = np.random.normal(900, 2, 3600)

# generate a star with a sine wave
flx_star = np.random.poisson(2700, 3600)
Fs = 3600
f = 5
sample = 3600
xx = np.arange(sample)
flx_var = ((np.sin(2 * np.pi * f * (xx / Fs)) * (2700 / 10)) + 2700) * (flx_star / 2700)
xx = np.arange(sample) * 3 / 60 / 60 / 24

# make a table of Gaussian sources
for ii in range(0, 3600):
    table = QTable()
    table['amplitude'] = [flx_var[ii], flx_star[ii]]
    table['x_mean'] = [x_var[ii], x_star[ii]]
    table['y_mean'] = [y_var[ii], y_star[ii]]
    table['x_stddev'] = [12, 12]
    table['y_stddev'] = [6, 6]
    table['theta'] = np.radians(np.array([0, 0]))

    image1 = make_gaussian_sources_image(shape, table)
    bkg = make_noise_image(shape, distribution='poisson', mean=5.)
    bias = make_noise_image(shape, distribution='gaussian', mean=95., stddev=0.5)

    image3 = image1 + bkg + bias

    if ii < 10:
        nme = '000' + str(ii)
    elif (ii >= 10) & (ii < 100):
        nme = '00' + str(ii)
    elif (ii >= 100) & (ii < 1000):
        nme = '0' + str(ii)
    else:
        nme = str(ii)

    min = ((ii * 3) // 60) % 60
    if min < 10:
        min_st = '0' + str(min)
    else:
        min_st = str(min)
    sec = (ii * 3) % 60
    if sec < 10:
        sec_st = '0' + str(sec)
    else:
        sec_st = str(sec)
    hr_st = '0' + str((ii * 3) // 3600)

    hdu = fits.PrimaryHDU()
    header = hdu.header

    header['DATE-OBS'] = '2022-08-21 ' + hr_st + ':' + min_st + ':' + sec_st + '.000000'
    header['HIERARCH EXPOSURE TIME'] = 3.0
    header['AIRMASS'] = 1.2
    header['FOCUS'] = 5100
    header['HIERARCH PRESSURE_MB'] = 900
    header['HIERARCH TEMPERATURE_C'] = 70
    header['HIERARCH HUMIDITY_PERCENT'] = 80
    header['HIERARCH DOME_HOUR'] = 10
    header['HIERARCH ZD_DEGREE'] = 90
    header['HIERARCH AZ_DEGREE'] = 90

    fits.writeto(dir + "test_image_" + nme + ".fits", image3, header=header, overwrite=True)

    if ii % 100 == 0:
        Utils.log("Working on next 100. " + str(3600 - ii - 1) + " images remain.", "info")

Utils.log("All fake images generated.", "info")
