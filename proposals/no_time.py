import pandas as pd

from libraries.utils import Utils
from astropy.io import fits
from photutils.centroids import centroid_sources, centroid_2dg
from photutils import aperture_photometry
from photutils.detection import find_peaks, IRAFStarFinder
from photutils import EllipticalAperture, EllipticalAnnulus, CircularAnnulus, CircularAperture
from photutils.psf import DAOGroup, BasicPSFPhotometry, extract_stars, EPSFBuilder, IterativelySubtractedPSFPhotometry
from photutils.background import MMMBackground, MADStdBackgroundRMS
from libraries.photometry import Photometry
from config import Configuration
from astropy.table import Table
import numpy as np
from astropy.stats import sigma_clipped_stats
nme = '20220609-TrES-1'
files = Utils.get_file_list('I:\\secondary_scope_images\\' + nme + '\\', ".FIT")

mags = np.zeros(len(files))
offs = np.zeros(len(files))
jd = np.zeros(len(files))
for idx, file in enumerate(files):
    img, header = fits.getdata('I:\\secondary_scope_images\\' + nme + '\\' + file, 0, header=True)

    # get the PSF peaks on the image
    peaks = find_peaks(img, threshold=Configuration.PSF_THRESHOLD)
    peaks = peaks[peaks['peak_value'] > 10000]

    star_list = pd.DataFrame(data=zip(peaks['x_peak'], peaks['y_peak']), columns=['x', 'y'])
    positions = star_list[['x', 'y']].apply(tuple, axis=1).tolist()

    # set up the apertures for the photometry
    aperture = CircularAperture(positions, 15)
    aperture_annulus = CircularAnnulus(positions, 17, 20)
    apers = [aperture, aperture_annulus]

    # run the photometry to get the data table
    phot_table = aperture_photometry(img, apers, method='exact')

    # subtract the sky background to get the stellar flux
    flx = np.array(phot_table['aperture_sum_0']) - ((phot_table['aperture_sum_1'] / aperture_annulus.area) * aperture.area)
    mag = 25 - 2.5 * np.log10(flx)

    if idx > 0:
        off = np.zeros(len(star_list))
        for ii in range(0, len(star_list)):
            dist = np.sqrt((star_list['x'][ii] - x0) ** 2 + (star_list['y'][ii] - y0) ** 2)
            minv = np.argmin(dist)
            off[ii] = mag0[minv] - mag[ii]
        _, mdn, stds = sigma_clipped_stats(off, sigma=2.5)
        offs[idx] = mdn
        mags[idx] = np.nanmin(mag)
    else:
        offs[idx] = 0
        mags[idx] = np.nanmin(mag)
        mag0 = mag
        x0 = star_list['x'].to_numpy()
        y0 = star_list['y'].to_numpy()

    jd[idx] = idx * 2

import matplotlib.pyplot as plt
plt.scatter(jd, mags)
plt.scatter(jd, mags+offs)
plt.show()
print('hold')