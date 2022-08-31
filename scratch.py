for idy, m in enumerate(mg):
    dt = {"xx": s[pos_x[idy]].to_numpy(), "yy": s[pos_y[idy]].to_numpy()}
    trnds = pd.DataFrame(dt)
    if (idy == 1) | (idy == 2) | (idy == 3) | (idy == 4):  ## | (idy == 5) | (idy == 6):

        for idx, file in enumerate(files):
            if idx == 0:
                continue
            else:
                tr = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + file, sep=' ', header=0)
                tr = tr.interpolate()
                tr['time_min'] = tr.apply(lambda x: (x.jd - tr.jd[0]) * 24. * 60., axis=1)
                # tr = tr[(tr.time_min < 13) | (tr.time_min > 31)].copy().reset_index(drop=True)
                # if s[m].corr(tr[m]) > 0.6:
                trnds['c' + str(idx)] = tr[m].to_numpy()
        # trnds['cc'] = s['mag_' + phot_type + '']
        # trnds = trnds.drop(columns=['xx', 'yy'])
        yhat = s[m].to_numpy()
        reg = LinearRegression().fit(trnds, yhat)
        trend = reg.predict(trnds)

        s[m] = s[m] - (trend - np.median(trend))
        s_bin = s.groupby('binned').agg({m: 'median', mg[1]: 'median'}).reset_index()
        er_bin = s.groupby('binned').agg({m: 'std', 'time_min': 'count'}).reset_index()

        plt.subplot(2, 1, 1)
        plt.errorbar(s_bin.binned, s_bin[m] - s_bin[m].median(), yerr=er_bin[m] / np.sqrt(er_bin.time_min),
                     fmt='none', c=clr[idy])
        plt.scatter(s_bin.binned, s_bin[m] - s_bin[m].median(), label=Configuration.WAVELENGTHS_TRANSMISSION[idy],
                    marker='o', c=clr[idy])
        plt.xlabel('Time from Observation Start [mins]')
        plt.ylabel('Relative Magnitude')
        plt.legend()
        # plt.ylim([0.03, -0.03])
        # plt.gca().invert_yaxis()

        plt.subplot(2, 1, 2)
        s['clr'] = s[m] - s[mg[1]]
        s_bin = s.groupby('binned').agg({'clr': 'median'}).reset_index()
        er_bin = s.groupby('binned').agg({'clr': 'std', 'time_min': 'count'}).reset_index()
        plt.errorbar(s_bin.binned, s_bin['clr'] - s_bin['clr'].median(),
                     yerr=er_bin['clr'] / np.sqrt(er_bin.time_min),
                     fmt='none', c=clr[idy])
        plt.scatter(s_bin.binned, s_bin['clr'] - s_bin['clr'].median(),
                    label=Configuration.WAVELENGTHS_TRANSMISSION[idy],
                    marker='o', c=clr[idy])
        plt.xlabel('Time from Observation Start [mins]')
        plt.ylabel('Relative Color [$\lambda$ - 763nm]')
        # plt.ylim([0.03, -0.03])
plt.gca().invert_yaxis()
plt.show()


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression
from config import Configuration
from libraries.utils import Utils


def custom_round(x, base=10):
    return int(base * round(float(x)/base))

if Configuration.STAR == 'XO-1':
    mg = ['mag_p1', 'mag_p2', 'mag_p3', 'mag_p4', 'mag_p5', 'mag_p6', 'mag_p7', 'mag_p8']
    pos_x = ['x_p1', 'x_p2', 'x_p3', 'x_p4', 'x_p5', 'x_p6', 'x_p7', 'x_p8']
    pos_y = ['y_p1', 'y_p2', 'y_p3', 'y_p4', 'y_p5', 'y_p6', 'y_p7', 'y_p8']

    clr = ['red', 'salmon', 'darkorange', 'gold', 'green', 'dodgerblue', 'blue', 'purple']

    s = pd.read_csv("F:\\ETSI_June2022\\2022-06-09\\T\\lc\\raw\\XO-1_transmission_00.lc", sep=' ',
                    header=0)
    s['time_min'] = s.apply(lambda x: (x.jd - s.jd[0]) * 24. * 60., axis=1)
    s['binned'] = s.apply(lambda x: custom_round(x.time_min), axis=1)
    s = s.interpolate()
    s = s[s.time_min < 400].copy().reset_index(drop=True)
    s17 = s.copy().reset_index()

    for idy, m in enumerate(mg):
        if idy < 5:
            dt = {"xx": s[pos_x[idy]].to_numpy(), "yy": s[pos_y[idy]].to_numpy()}
            trnds = pd.DataFrame(dt)

            tr = pd.read_csv("F:\\ETSI_June2022\\2022-06-17\\XO-1\\lc\\raw\\XO-1_transmission_01.lc",
                             sep=' ', header=0)
            tr['time_min'] = tr.apply(lambda x: (x.jd - tr.jd[0]) * 24. * 60., axis=1)
            tr = tr[tr.time_min < 400].copy().reset_index(drop=True)
            tr = tr.interpolate()
            trnds['c'] = tr[m].to_numpy()

            # trnds = trnds.drop(columns=['xx', 'yy'])
            yhat = s[m].to_numpy()
            reg = LinearRegression().fit(trnds, yhat)
            trend = reg.predict(trnds)

            s17[m] = s[m] - (trend - np.median(trend))

    # plt.subplot(2, 1, 1)
    for idy, m in enumerate(mg):
        if idy < 5:
            plt.subplot(2, 1, 1)
            s17_bin = s17.groupby('binned').agg({m: 'median'}).reset_index()
            er17_bin = s17.groupby('binned').agg({m: 'std', 'time_min': 'count'}).reset_index()

            plt.scatter(s17_bin['binned'], s17_bin[m] - s17_bin[m].median()+.01, c=clr[idy],
                        label=Configuration.WAVELENGTHS_TRANSMISSION[idy], marker='o')
            plt.errorbar(s17_bin['binned'], s17_bin[m] - s17_bin[m].median()+.01, yerr=er17_bin[m]/er17_bin['time_min'],
                         c=clr[idy], marker='o', fmt='none')

            plt.xlabel('Time from initial exposure [min]')
            plt.ylabel('Relative Magnitude')
            plt.legend()
            plt.ylim([0.04, -0.02])

    for idy, m in enumerate(mg):
        if idy < 5:
            plt.subplot(2, 1, 2)
            s17['clr'] = (s17[m]-s17[m].median()) - (s17[mg[2]]-s17[mg[2]].median())
            s_bin = s17.groupby('binned').agg({'clr': 'median'}).reset_index()
            er_bin = s17.groupby('binned').agg({'clr': 'std', 'time_min': 'count'}).reset_index()
            plt.errorbar(s_bin.binned, s_bin['clr'] - s_bin['clr'].median(),
                         yerr=er_bin['clr'] / np.sqrt(er_bin.time_min),
                         fmt='none', c=clr[idy])
            plt.scatter(s_bin.binned, s_bin['clr'] - s_bin['clr'].median(),
                        label=Configuration.WAVELENGTHS_TRANSMISSION[idy],
                        marker='.', c=clr[idy])
            plt.xlabel('Time from Observation Start [mins]')
            plt.ylabel('Relative Color [$\lambda$ - 660nm, mag]')
            plt.ylim([0.01, -0.01])
            plt.legend()
plt.show()

if Configuration.STAR == 'NGC6886':
    s_bin = s[mg].agg({'median', 'std', 'count'})

    plt.subplot(2, 1, 1)
    plt.errorbar(Configuration.WAVELENGTHS_TRANSMISSION_NUMS, s_bin.loc['median'].to_numpy()-s_bin.loc['median'][1],
                 yerr=s_bin.loc['std'].to_numpy()/s_bin.loc['count'].to_numpy(), fmt='none', color='k')
    plt.scatter(Configuration.WAVELENGTHS_TRANSMISSION_NUMS,
                s_bin.loc['median'].to_numpy()-s_bin.loc['median'][1], color='k')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Relative Color [$\lambda$ - 763nm, mag]')
    plt.legend()
    plt.gca().invert_yaxis()

    for idy, m in enumerate(mg):
        # if (idy == 1) | (idy == 2) | (idy == 3) | (idy == 4) | (idy == 5):
        plt.subplot(2, 1, 2)
        s['clr'] = s[m] - s[mg[1]]
        s_bin = s.groupby('binned').agg({'clr': 'median'}).reset_index()
        er_bin = s.groupby('binned').agg({'clr': 'std', 'time_min': 'count'}).reset_index()
        plt.errorbar(s_bin.binned, s_bin['clr'] - s_bin['clr'].median(),
                     yerr=er_bin['clr'] / np.sqrt(er_bin.time_min),
                     fmt='none', c=clr[idy])
        plt.scatter(s_bin.binned, s_bin['clr'] - s_bin['clr'].median(),
                    label=Configuration.WAVELENGTHS_TRANSMISSION[idy],
                    marker='.', c=clr[idy])
        plt.xlabel('Time from Observation Start [mins]')
        plt.ylabel('Relative Color [$\lambda$ - 763nm, mag]')
        plt.ylim([0.004, -0.004])
        plt.legend()
    plt.show()


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


files = Utils.get_file_list(Configuration.LIGHTCURVE_DIRECTORY, '.lc')

nfiles = len(files)
ntrnds = len(files) - 1

s = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + files[0], sep=' ', header=0)
s['time_min'] = np.arange(98) * 1.34  # s.apply(lambda x: (x.jd - s.jd[0]) * 24. * 60., axis=1)
s['binned'] = s.apply(lambda x: custom_round(x.time_min), axis=1)
s = s.interpolate()
# s = s[(s.time_min < 13) | (s.time_min > 31)].copy().reset_index(drop=True)

mg = ['mag_p1', 'mag_p2', 'mag_p3', 'mag_p4', 'mag_p5', 'mag_p6', 'mag_p7', 'mag_p8']
pos_x = ['x_p1', 'x_p2', 'x_p3', 'x_p4', 'x_p5', 'x_p6', 'x_p7', 'x_p8']
pos_y = ['y_p1', 'y_p2', 'y_p3', 'y_p4', 'y_p5', 'y_p6', 'y_p7', 'y_p8']

clr = ['red', 'salmon', 'darkorange', 'gold', 'green', 'dodgerblue', 'blue', 'purple']

for idy, m in enumerate(mg):
    dt = {"xx": s[pos_x[idy]].to_numpy(), "yy": s[pos_y[idy]].to_numpy()}
    trnds = pd.DataFrame(dt)

    for idx, file in enumerate(files):
        if idx == 0:
            continue
        else:
            tr = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + file, sep=' ', header=0)
            tr = tr.interpolate()
            # tr['time_min'] = tr.apply(lambda x: (x.jd - tr.jd[0]) * 24. * 60., axis=1)
            # tr = tr[(tr.time_min < 13) | (tr.time_min > 31)].copy().reset_index(drop=True)
            trnds['c' + str(idx)] = tr[m].to_numpy()
    trnds = trnds.drop(columns=['xx', 'yy'])
    yhat = s[m].to_numpy()
    reg = LinearRegression().fit(trnds, yhat)
    trend = reg.predict(trnds)

    s[m] = s[m] - (trend - np.median(trend))


def custom_round(x, base=10):
    return int(base * round(float(x)/base))

date = '2022-06-17'
s22 = pd.read_csv("F:\\ETSI_June2022\\" + date + "\\XO-1\\lc\\raw\\og_bkg\\XO-1_transmission_00.lc", sep=' ', header=0)
c22 = pd.read_csv("F:\\ETSI_June2022\\" + date + "\\XO-1\\lc\\raw\\og_bkg\\XO-1_transmission_01.lc", sep=' ', header=0)

date = '2022-06-17'
s17 = pd.read_csv("F:\\ETSI_June2022\\" + date + "\\XO-1\\lc\\raw\\XO-1_transmission_00.lc", sep=' ', header=0)
c17 = pd.read_csv("F:\\ETSI_June2022\\" + date + "\\XO-1\\lc\\raw\\XO-1_transmission_01.lc", sep=' ', header=0)

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