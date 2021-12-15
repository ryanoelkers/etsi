from config import Configuration
from astropy.io import fits
import numpy as np
import matplotlib
import pandas as pd
from sklearn import metrics
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
from scipy.spatial.distance import cdist
from astropy.stats import sigma_clipped_stats
from sklearn.cluster import KMeans
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn import linear_model


images = fits.getdata(Configuration.TEST_DIRECTORY + 'HD213786_16bit_drift_0p5s_1.fits', 0)

mx_img = np.amax(images[0:27, :, :], axis=0) / np.max(images[0:27, :, :])
mx_hdu = fits.PrimaryHDU(data=mx_img)
mx_hdu.writeto(Configuration.ANALYSIS_DIRECTORY + 'max_27_image.fits', overwrite=True)

mx_img = np.amax(images[121:144, :, :], axis=0) / np.max(images[121:144, :, :])
mx_hdu = fits.PrimaryHDU(data=mx_img)
mx_hdu.writeto(Configuration.ANALYSIS_DIRECTORY + 'max_121_image.fits', overwrite=True)

# get the max image
mx_img = np.amax(images, axis=0) / np.max(images)
mx_hdu = fits.PrimaryHDU(data=mx_img)
mx_hdu.writeto(Configuration.ANALYSIS_DIRECTORY + 'max_image.fits', overwrite=True)

# get the min image
mn_img = np.amin(images, axis=0)
mn_hdu = fits.PrimaryHDU(data=mn_img)
mn_hdu.writeto(Configuration.ANALYSIS_DIRECTORY + 'min_image.fits', overwrite=True)

# get the median image
md_img = np.median(images, axis=0)
md_hdu = fits.PrimaryHDU(data=md_img)
md_hdu.writeto(Configuration.ANALYSIS_DIRECTORY + 'median_image.fits', overwrite=True)

for image in images:

    # get the background on the image
    bkg = np.median(image)

    # get a background cutoff
    mn_bkg, md_bkg, sg_bkg = sigma_clipped_stats(image, sigma_upper=2)
    cut_bkg = mn_bkg + 50 * sg_bkg
    X = np.argwhere(image > cut_bkg)

    kmeans = KMeans(n_clusters=4, random_state=0).fit(X)
    positions = pd.DataFrame(kmeans.cluster_centers_, columns=['y', 'x'])
    pixels = pd.DataFrame(X, columns=['y', 'x'])
    pixels['labs'] = kmeans.labels_

    # get the pixel count for each star
    stars = pixels.groupby('labs').agg({'x': 'count'})
    stars = stars.rename(columns={'x': 'pixel_count'})
    stars['x'] = kmeans.cluster_centers_[:, 1]
    stars['y'] = kmeans.cluster_centers_[:, 0]

    stp = 0
    cutoff = mn_bkg + 5 * sg_bkg
    stpixs = stars[['y', 'x']].astype(int).to_numpy()

airmass = np.array([1.05, 1.08])
s = 0.1 * (210. ** (-2./3.)) * (airmass ** (7./4.)) * np.exp(-2070./8000.) * ((2.*0.25) ** (-1./2.))
lc = pd.read_csv('/home/oelkerrj/Research/ETSI/testing/code/data.txt', header=0)
reg = linear_model.LinearRegression()

lc['dist'] = lc.apply(lambda x: np.sqrt((x.x - lc.x[0])**2 + (x.y - lc.y[0])**2), axis=1)

X = lc[['trd']].to_numpy().reshape(-1, 1)
y = lc['mag'].to_numpy()
model = reg.fit(X, y)

reg_mod = model.coef_[0] * lc.trd.to_numpy() + model.intercept_
lc['cln'] = lc.mag-reg_mod

plt.figure(figsize=(7, 11))
plt.subplot(3, 1, 1)
plt.scatter(lc.jd, lc.mag)
plt.ylim([0.03, -0.03])
plt.ylabel('Relative Magnitude Primary [mag]')
plt.xlabel('JD - 2459465 [d]')
# plt.title('RMS ~ 0.0054 mag')

plt.subplot(3, 1, 2)
plt.scatter(lc.jd, lc.trd)
plt.ylim([0.03, -0.03])
plt.ylabel('Relative Magnitude Secondary [mag]')
plt.xlabel('JD - 2459465 [d]')
# plt.title('RMS ~ 0.0085 mag')

plt.subplot(3, 1, 3)
plt.scatter(lc.jd, lc.mag - lc.trd)
plt.ylim([0.03, -0.03])
plt.ylabel('Primary - Secondary [mag]')
plt.xlabel('JD - 2459465 [d]')
# plt.title('RMS ~ 0.0076 mag')

plt.savefig('/home/oelkerrj/Research/ETSI/testing/code/lc_detrend.png', bbox='tight')

plt.figure(figsize=(7, 11))
plt.subplot(3, 1, 1)
plt.scatter(lc.jd, lc.mag, label='primary')
plt.plot(lc.jd, reg_mod, label='model', c='red')
plt.ylim([0.02, -0.02])
plt.ylabel('Relative Magnitude [mag]')
plt.xlabel('JD - 2459465 [d]')
# plt.title('RMS ~ 0.0054 mag')

plt.subplot(3, 1, 2)
plt.scatter(lc.jd, lc.cln)
plt.ylim([0.02, -0.02])
plt.ylabel('Relative Magnitude [mag]')
plt.xlabel('JD - 2459465 [d]')
# plt.title('RMS ~ 0.0047 mag')


plt.subplot(3, 1, 3)
plt.scatter(lc[(lc.cln < 0.01) & (lc.cln > -0.01)].x,
            lc[(lc.cln < 0.01) & (lc.cln > -0.01)].y,
            c=lc[(lc.cln < 0.01) & (lc.cln > -0.01)].cln)
plt.colorbar()
plt.savefig('/home/oelkerrj/Research/ETSI/testing/code/pos_lc_detrend.png', bbox='tight')

images = fits.getdata(Configuration.TEST_DIRECTORY + 'HD213786_16bit_drift_0p5s_1.fits', 0)

for image in images:

    # get the background on the image
    bkg = np.median(image)

    # get a background cutoff
    mn_bkg, md_bkg, sg_bkg = sigma_clipped_stats(image, sigma_upper=2)
    cut_bkg = mn_bkg + 50 * sg_bkg
    X = np.argwhere(image > cut_bkg)

    kmeans = KMeans(n_clusters=4, random_state=0).fit(X)
    positions = pd.DataFrame(kmeans.cluster_centers_, columns=['y', 'x'])
    pixels = pd.DataFrame(X, columns=['y', 'x'])
    pixels['labs'] = kmeans.labels_

    # get the pixel count for each star
    stars = pixels.groupby('labs').agg({'x': 'count'})
    stars = stars.rename(columns={'x': 'pixel_count'})
    stars['x'] = kmeans.cluster_centers_[:, 1]
    stars['y'] = kmeans.cluster_centers_[:, 0]

    stp = 0
    cutoff = mn_bkg + 5 * sg_bkg
    stpixs = stars[['y', 'x']].astype(int).to_numpy()

    plt.scatter(pixels.x, pixels.y, c='r', marker='.')
    plt.scatter(positions.x, positions.y, c='b', marker='x')
    plt.xlim([0, 2048])
    plt.ylim([0, 2048])
    plt.xlabel('X Pixel')
    plt.ylabel('Y Pixel')
    plt.show()



    for star in star_idx:
        plt.scatter([star[0]], [star[1]], marker='x', c='.')

    plt.show()
    print('hold')
