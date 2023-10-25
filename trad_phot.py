import numpy as np
import pandas as pd
import matplotlib
from sklearn.linear_model import LinearRegression
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from config import Configuration
import os

mg = ['aper_1', 'aper_2', 'aper_3', 'aper_4', 'aper_5', 'aper_6', 'aper_7', 'aper_8']
cmg = ['aper_c1', 'aper_c2', 'aper_c3', 'aper_c4', 'aper_c5', 'aper_c6', 'aper_c7', 'aper_c8']
mg_e = ['aper_e1', 'aper_e2', 'aper_e3', 'aper_e4', 'aper_e5', 'aper_e6', 'aper_e7', 'aper_e8']
bkg = ['aper_bkg1', 'aper_bkg2', 'aper_bkg3', 'aper_bkg4', 'aper_bkg5', 'aper_bkg6', 'aper_bkg7', 'aper_bkg8']
clrt = ['maroon', 'crimson', 'goldenrod', 'yellowgreen', 'blue', 'purple']
clrr = ['red', 'orange', 'gold', 'green', 'dodgerblue']
# get the light curve in transmission and reflection
s = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + Configuration.STAR + '_transmission_00_none.lc',
                sep=' ', header=0)
c = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + Configuration.STAR + '_transmission_01_none.lc',
                sep=' ', header=0)
r = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + Configuration.STAR + '_reflection_00_none.lc',
                sep=' ', header=0)
t = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + Configuration.STAR + '_reflection_01_none.lc',
                sep=' ', header=0)

# subtract the background
for idx, m in enumerate(mg):

    s[mg_e[idx]] = s.apply(lambda x: x[mg_e[idx]] / (x[m] - x[bkg[idx]]), axis=1)
    s[m] = s.apply(lambda x: 25 - 2.5 * np.log10((x[m] - x[bkg[idx]]) / x['exp_time']), axis=1)
    s['tr'] = c.apply(lambda x: 25 - 2.5 * np.log10((x[m] - x[bkg[idx]]) / x['exp_time']), axis=1)

    xx = s[['tr']]
    yy = s[m]
    reg = LinearRegression().fit(xx, yy).predict(xx)

    s[cmg[idx]] = s[m] # - reg

    if idx < 7:
        r[mg_e[idx]] = r.apply(lambda x: x[mg_e[idx]] / (x[m] - x[bkg[idx]]), axis=1)
        r[m] = r.apply(lambda x: 25 - 2.5 * np.log10((x[m] - x[bkg[idx]]) / x['exp_time']), axis=1)
        r['tr'] = t.apply(lambda x: 25 - 2.5 * np.log10((x[m] - x[bkg[idx]]) / x['exp_time']), axis=1)

        xx = r[['tr']]
        yy = r[m]
        reg = LinearRegression().fit(xx, yy).predict(xx)

        r[cmg[idx]] = r[m] # - reg

plt.figure(figsize=(9, 6))

for idx, m in enumerate(cmg):
    if idx < 6:
        # plot the light curve
        plt.errorbar(s.airmass, s[m] - s[cmg[2]] - np.median(s[m] - s[cmg[2]]),
                     yerr=np.sqrt(s[mg_e[idx]] ** 2 + s[mg_e[2]] ** 2), c=clrt[idx], fmt='none')

        plt.plot(s.airmass, s[m] - s[cmg[2]] - np.median(s[m] - s[cmg[2]]),
                 c=clrt[idx], label=Configuration.WAVELENGTHS_TRANSMISSION[idx], marker='.')

    if idx < -1:
        # plot the light curve
        plt.errorbar((r.jd - r.jd.min()) * 24, r[m] - s[cmg[2]] + (0.02 - (idx * 2 + 1) * 0.005),
                     yerr=np.sqrt(r[mg_e[idx]] ** 2 + s[mg_e[2]] ** 2), c=clrr[idx], fmt='none')

        plt.plot((r.jd - r.jd.min()) * 24, r[m] - s[cmg[2]] + (0.02 - (idx * 2 + 1) * 0.005),
                 c=clrr[idx], label=Configuration.WAVELENGTHS_REFLECTION[idx], marker='.')

plt.xlabel('Time Since First Exposure [hours]', fontsize=15)
plt.ylabel('Relative Color [mag, $\lambda$ - 660nm]', fontsize=15)
plt.legend(loc='upper left', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().invert_yaxis()
# plt.savefig("C:\\Users\\barristan\\Desktop\\BL_LAC.png", bbox_inches='tight', pad_inches=0.1)
plt.show()

print('hold')