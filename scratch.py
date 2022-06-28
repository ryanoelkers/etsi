""" This set of functions is primarily used for photometery."""
from config import Configuration
import pandas as pd
import numpy as np
import warnings
from sklearn.linear_model import LinearRegression
from scipy.signal import savgol_filter
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import r2_score
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import logging
logging.getLogger('matplotlib.font_manager').disabled = True


class Lightcurves:

    @staticmethod
    def detrend_lightcurves(lightcurve_directory, ):
        """ This function will detrended the target star light curve using a comparison star or a synthetic trend if
        no comparison star is available for the detrending.

        :parameter lightcurve_directory - The directory with the light curves

        :return Nothing is returned, but the detrended data are printed to the cleaned directory
        """

        def custom_round(x, base=10):
            return int(base * round(float(x) / base))
        # get the target star light curve, should always be star 00
        lc = pd.read_csv("C:\\Users\\barristan\\Documents\\ETSI\\data\\HAT-P-44\\lc\\raw\\HAT-P-44_00.lc", sep=' ')
        tr = pd.read_csv("C:\\Users\\barristan\\Documents\\ETSI\\data\\HAT-P-44\\lc\\raw\\HAT-P-44_01.lc", sep=' ')
        lc = lc.interpolate()
        lc = lc[(lc.j >= 20) & (lc.j <= 110)].copy().reset_index(drop=True)
        lc['j'] = lc['j'] - 20
        lc['binned'] = lc.apply(lambda x: custom_round(x.j), axis=1)

        mg = ['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8']
        tg = ['tm1', 'tm2', 'tm3', 'tm4', 'tm5', 'tm6', 'tm7', 'tm8']
        ms = ['s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8']
        ts = ['ts1', 'ts2', 'ts3', 'ts4', 'ts5', 'ts6', 'ts7', 'ts8']
        pg = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8']

        st_mod = lc[mg].median(axis=1)
        err = np.zeros(8)
        for idx, m in enumerate(mg):

            lcm = savgol_filter(lc[m], 1501, 3, mode='nearest')
            str = savgol_filter(st_mod, 1501, 1, mode='nearest')

            # trd = st_mod.rolling(1500, center=True, min_periods=1).mean()
            # trd = st_mod
            y = lc[m].to_numpy()
            x = np.array(str).reshape(-1, 1)

            # kernel = DotProduct() + WhiteKernel()
            # gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(x, y)
            polynomial_features = PolynomialFeatures(degree=2)
            x_poly = polynomial_features.fit_transform(x)

            model = LinearRegression()
            model.fit(x_poly, y)

            yhat = model.predict(x_poly)
            # reg = LinearRegression().fit(x, y)
            # yhat = reg.predict(x)

            lc[pg[idx]] = y - (yhat - np.median(yhat))
            # lc[pg[idx]] = y - (str - np.median(str))
            if idx == 0:
                yhat2 = savgol_filter(lc[pg[idx]], 5001, 1, mode='nearest')
                lc[pg[idx]] = lc[pg[idx]] - (yhat2 - np.median(yhat2))
            #plt.scatter(lc.j, lc[m])
            # plt.scatter(lc.j, trd - np.median(trd) + np.median(lc[m]), marker='.')
            #plt.scatter(lc.j, yhat - np.median(yhat) + np.median(lc[m]), marker='.')
            #plt.scatter(lc.j, lcm - np.median(lcm) + np.median(lc[m]), marker='.')
            # plt.scatter(lc.j, y - (trd - np.median(trd)), marker='.')
            # plt.show()
            flx1 = lc.apply(lambda x: 10 ** ((x[m] - 25)/(-2.5)), axis=1) * 0.61
            err[idx] = np.mean(np.sqrt(flx1) / (flx1 - (np.pi * 35 * 15 * 100 * 0.61))) / np.sqrt(3000)

        lc_bin = lc.groupby('binned').agg({'p1': 'mean', 'p2': 'mean',
                                           'p3': 'mean', 'p4': 'mean',
                                           'p5': 'mean', 'p6': 'mean',
                                           'p7': 'mean', 'p8': 'mean'}).reset_index()
        lc_std = lc.groupby('binned').agg({'p1': 'std', 'p2': 'std',
                                           'p3': 'std', 'p4': 'std',
                                           'p5': 'std', 'p6': 'std',
                                           'p7': 'std', 'p8': 'std'}).reset_index()
        lc_cnt = lc.groupby('binned').agg({'p1': 'count', 'p2': 'count',
                                           'p3': 'count', 'p4': 'count',
                                           'p5': 'count', 'p6': 'count',
                                           'p7': 'count', 'p8': 'count'}).reset_index()

        cc = 'p3'
        clr = ['maroon', 'red', 'salmon', 'darkorange', 'gold', 'darkgreen', 'blue', 'purple']
        plt.figure(figsize=(12, 8), dpi=120)
        for idx, m in enumerate(pg):
            # trd = lc['p3'].rolling(50, center=True, min_periods=1).mean()
            plt.scatter(lc_bin.binned, lc_bin[m] - lc_bin[cc] - np.median(lc_bin[m] - lc_bin[cc]) + (-1*idx*0.005),
                        label=Configuration.WAVELENGTHS_TRANSMISSION[idx], c=clr[idx])
            c = np.polyfit(lc_bin.binned, lc_bin[m] - lc_bin[cc] - np.median(lc_bin[m] - lc_bin[cc]), 1)
            yhat = c[0] * lc_bin.binned + c[1]
            # print(r2_score(lc_bin[m] - lc_bin[cc], yhat))
            plt.plot(lc_bin.binned, 0 * lc_bin.binned + (-1*idx * 0.005), c=clr[idx], linewidth=2)
            # plt.plot(lc_bin.binned, yhat + (idx * 0.003), c=clr[idx])
            plt.errorbar(lc_bin.binned, lc_bin[m] - lc_bin[cc] - np.median(lc_bin[m] - lc_bin[cc]) + (-1*idx*0.005),
                         yerr=lc_std[m] / (lc_std[m] / lc_bin[m].std()), fmt='none', c=clr[idx])  # lc_cnt[m]
            print(np.std(lc_bin[m] - lc_bin[cc]), lc_std[m].mean() / lc_cnt[m].mean(), err[idx])
        plt.ylabel('Relative Color (- offset) [$\lambda$ - 660nm] [mag]', size=17)
        plt.xlabel('Time Since First Observation [mins]', size=17)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.ylim([-0.04, 0.005])
        plt.legend(loc='upper right', fontsize=15)
        plt.savefig(Configuration.ANALYSIS_DIRECTORY + 'HD94883.png', bbox_inches='tight', pad_inches=0)
        #plt.show()
        print('hold')

@staticmethod
def mk_flat(image_directory):
    """ This function will make the master flat frame using the provided image list.
    :parameter image_directory - a directory where the images reside for combination

    :return - The bias frame is returned and written to the master directory
    """

    if os.path.isfile(Configuration.CALIBRATION_DATE_DIRECTORY + Configuration.DATE + "_flat.fits") == 0:

        # read in the bias file
        bias = Preprocessing.mk_bias(Configuration.BIAS_DIRECTORY)

        # get the image list
        image_list = Utils.get_file_list(image_directory, '.fits')

        # determine the number of loops we need to move through for each image
        nfiles = len(image_list)
        nbulk = 20

        # get the integer and remainder for the combination
        full_bulk = nfiles // nbulk
        part_bulk = nfiles % nbulk
        if part_bulk > 0:
            hold_bulk = full_bulk + 1
        else:
            hold_bulk = full_bulk

        # here is the 'holder'
        hold_data = np.array(shape=(Configuration.AXS, Configuration.AXS))
        # np.ndarray(shape=(hold_bulk, Configuration.AXS, Configuration.AXS))

        # update the log
        Utils.log("Generating a master flat frame from multiple files in bulks of " + str(nbulk) +
                  " images. There are " + str(nfiles) + " images to combine, which means there should be " +
                  str(hold_bulk) + " mini-files to combine.", "info")

        for kk in range(0, hold_bulk):

            # loop through the images in sets of nbulk
            if kk < full_bulk:
                # generate the image holder
                block_hold = np.ndarray(shape=(nbulk, Configuration.AXS, Configuration.AXS))

                # generate the max index
                mx_index = nbulk
            else:
                # generate the image holder
                block_hold = np.ndarray(shape=(part_bulk, Configuration.AXS, Configuration.AXS))

                # generate the max index
                mx_index = part_bulk

            # make the starting index
            loop_start = kk * nbulk
            idx_cnt = 0

            Utils.log("Making mini-flat file " + str(kk) + ".", "info")

            # now loop through the images
            for jj in range(loop_start, mx_index + loop_start):
                # read in the flat file
                flat_tmp = fits.getdata(Configuration.FLAT_DIRECTORY + image_list[jj], 0)
                if np.ndim(flat_tmp) > 2:
                    flat_tmp = flat_tmp[0]

                # read in the image directly into the block_hold
                block_hold[idx_cnt] = flat_tmp - bias

                # increase the iteration
                idx_cnt += 1

            # median the data into a single file
            hold_data[kk] = np.median(block_hold, axis=0)

        # median the mini-images into one large image
        flat_image = np.median(hold_data, axis=0)
        nflat_image = flat_image / np.median(flat_image[512:1536, 512:1536])

        # pull the header information from the first file of the set
        flat_header = fits.getheader(Configuration.FLAT_DIRECTORY + image_list[0])
        flat_header['comb_typ'] = 'median'
        flat_header['norm_pix'] = '512to1536'
        flat_header['num_comb'] = nfiles

        # write the image out to the master directory
        fits.writeto(Configuration.CALIBRATION_DATE_DIRECTORY + Configuration.DATE + "_flat.fits",
                     nflat_image, flat_header, overwrite=True)
    else:
        nflat_image = fits.getdata(Configuration.CALIBRATION_DATE_DIRECTORY + Configuration.DATE + "_flat.fits", 0)

    return nflat_image