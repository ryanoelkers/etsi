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
        lc = pd.read_csv("C:\\Users\\barristan\\Documents\\ETSI\\data\\HAT-P-44\\lc\\raw\\sky\\HAT-P-44_00.lc", sep=' ')
        tr = pd.read_csv("C:\\Users\\barristan\\Documents\\ETSI\\data\\HAT-P-44\\lc\\raw\\sky\\HAT-P-44_01.lc", sep=' ')
        lc = lc.interpolate()

        mg = ['mag_p1', 'mag_p2', 'mag_p3', 'mag_p4', 'mag_p5', 'mag_p6', 'mag_p7', 'mag_p8']

        for idx, m in enumerate(mg):

            trd = tr[m].rolling(151, center=True, min_periods=1).mean()
            y = lc[m].to_numpy()
            x = np.array(tr[m]).reshape(-1, 1)

            reg = LinearRegression().fit(x, y)
            yhat = reg.predict(x)

            # lc[pg[idx]] = y - (yhat - np.median(yhat))

            plt.scatter(lc.jd, lc[m] - (yhat - np.median(yhat) + np.median(lc[m])), marker='.')
            plt.show()

        print('hold')