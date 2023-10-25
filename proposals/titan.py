import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, medfilt
from sklearn.linear_model import LinearRegression
from config import Configuration
from scipy import interpolate
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import r2_score
from astropy.stats import sigma_clipped_stats


# get the names of the relevant files
mg = ['aper_1', 'aper_2', 'aper_3', 'aper_4', 'aper_5', 'aper_6', 'aper_7', 'aper_8']
psf = ['psf_1', 'psf_2', 'psf_3', 'psf_4', 'psf_5', 'psf_6', 'psf_7', 'psf_8']

mxs = ['aper_max1', 'aper_max2', 'aper_max3', 'aper_max4', 'aper_max5', 'aper_max6', 'aper_max7', 'aper_max8']

nmg = ['aper_n1', 'aper_n2', 'aper_n3', 'aper_n4', 'aper_n5', 'aper_n6', 'aper_n7', 'aper_n8']
cmg = ['aper_c1', 'aper_c2', 'aper_c3', 'aper_c4', 'aper_c5', 'aper_c6', 'aper_c7', 'aper_c8']
clrmg = ['aper_clr1', 'aper_clr2', 'aper_clr3', 'aper_clr4', 'aper_clr5', 'aper_clr6', 'aper_clr7', 'aper_clr8']
bkg = ['aper_bkg1', 'aper_bkg2', 'aper_bkg3', 'aper_bkg4', 'aper_bkg5', 'aper_bkg6', 'aper_bkg7', 'aper_bkg8']
skval = ['aper_sky1', 'aper_sky2', 'aper_sky3', 'aper_sky4', 'aper_sky5', 'aper_sky6', 'aper_sky7', 'aper_sky8']

xx = ['x_a1', 'x_a2', 'x_a3', 'x_a4', 'x_a5', 'x_a6', 'x_a7', 'x_a8']
yy = ['y_a1', 'y_a2', 'y_a3', 'y_a4', 'y_a5', 'y_a6', 'y_a7', 'y_a8']

bias_level = 98. * (np.pi * Configuration.ELIP_APER_B * Configuration.ELIP_APER_A) * \
             Configuration.BIN_TIME / Configuration.EXPOSURE_TIME
bias_noise = 3.5 * (np.pi * Configuration.ELIP_APER_B * Configuration.ELIP_APER_A) * \
             Configuration.BIN_TIME / Configuration.EXPOSURE_TIME
npix = np.pi * Configuration.ELIP_APER_B * Configuration.ELIP_APER_A

clr = ['maroon', 'red', 'darkorange', 'orange', 'black', 'gold', 'yellow',
       'darkgreen', 'lime', 'teal', 'dodgerblue', 'blue', 'indigo',
       'purple', 'fuchsia']

# get the light curve
s_tr = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + Configuration.STAR + '_transmission_00_none.lc',
                   sep=' ', header=0)
s_tr['time_min'] = s_tr.apply(lambda x: (x['jd'] - s_tr['jd'].min()) * 24 * 60, axis=1)
s_rf = pd.read_csv(Configuration.LIGHTCURVE_DIRECTORY + Configuration.STAR + '_reflection_00_none.lc',
                   sep=' ', header=0)
s_rf['time_min'] = s_tr.apply(lambda x: (x['jd'] - s_rf['jd'].min()) * 24 * 60, axis=1)
# add minutes from start
tr_st = 2459769.88552
tr_ed = 2459769.89393

spectra = pd.read_csv(Configuration.ANALYSIS_DIRECTORY + Configuration.STAR + '\\1995low.tab', sep=r"\s+")
wvs = np.sort(np.concatenate([np.array(Configuration.WAVELENGTHS_TRANSMISSION_NUMS),
                              np.array(Configuration.WAVELENGTHS_REFLECTION_NUMS)]))
spect = np.interp(wvs, spectra.wve_vc, spectra.titan)
spect_nrm = spect / spect[10]

wvs_up = [900, 847, 792.5, 697, 663.5, 600.5,   574.5, 544,   523, 503, 485.5, 468,   456.5, 440,   430]
wvs_dw = [975, 900, 897,   729.5, 697,   630,   600.5, 574.5, 544, 523, 503,   485.5, 468,   456.5, 440]

depths = np.zeros(15)
ers = np.zeros(15)

spectra2 = pd.read_csv(Configuration.ANALYSIS_DIRECTORY + Configuration.STAR + '\\titan_spectra_2.txt', sep=r"\s+")
spect2 = np.interp(wvs, spectra2.wve/10, spectra2.flx)

r = np.mean([8.925652, 8.938952])
rp = 2560.0
delta = np.mean([9.56266, 9.56298])

spectra3 = pd.read_csv(Configuration.ANALYSIS_DIRECTORY + Configuration.STAR + '\\spectrum3.txt', sep=r"\s+")
spectra3['flx'] = spectra3.apply(lambda x:
                                 (((x.titan / 1000.) * (rp / delta)**2) / (r**2)) * x.sun, axis=1)
spect3 = np.interp(wvs, spectra3.wve/10, spectra3.flx)

for idx, wv in enumerate(wvs):
    depths[idx] = spectra[(spectra['wve_vc'] > wvs_up[idx]) & (spectra['wve_vc'] < wvs_dw[idx])]['titan'].mean()
    ers[idx] = spectra[(spectra['wve_vc'] > wvs_up[idx]) & (spectra['wve_vc'] < wvs_dw[idx])]['titan'].std()

# subtract the background
for idx, b in enumerate(bkg):
    s_tr[mg[idx]] = s_tr.apply(lambda x: (x[mg[idx]] - x[b]) / x['exp_time'], axis=1)

    if idx < 7:
        s_rf[mg[idx]] = s_rf.apply(lambda x: (x[mg[idx]] - x[b]) / x['exp_time'], axis=1)

for idx, b in enumerate(bkg):
    knot_numbers = 5
    x_new = np.linspace(0, 1, knot_numbers + 2)[1:-1]
    q_knots = np.quantile(s_tr.time_min.to_numpy(), x_new)
    #
    tt, c, k = interpolate.splrep(s_tr[(s_tr.jd < tr_st) | (s_tr.jd > tr_ed)].time_min.to_numpy(),
                                  s_tr[(s_tr.jd < tr_st) | (s_tr.jd > tr_ed)][mg[idx]].to_numpy(),
                                  t=q_knots, s=1)

    yfit = interpolate.BSpline(tt, c, k)(s_tr.time_min.to_numpy())
    s_tr[cmg[idx]] = s_tr[mg[idx]] / yfit

    # plt.subplot(2, 1, 1)
    # plt.scatter(s_tr.time_min, s_tr[mg[idx]], c='k')
    # plt.plot(s_tr.time_min, yfit, c='b')

    # plt.subplot(2, 1, 2)
    # plt.plot(s_tr.time_min, s_tr[cmg[idx]], c=clr[idx * 2],
    #         label='Transmission ' + Configuration.WAVELENGTHS_TRANSMISSION[idx])
    # plt.legend()
    # plt.show()

    if idx < 7:
        knot_numbers = 5
        x_new = np.linspace(0, 1, knot_numbers + 2)[1:-1]
        q_knots = np.quantile(s_rf.time_min.to_numpy(), x_new)
        tt, c, k = interpolate.splrep(s_rf[(s_rf.jd < tr_st) | (s_rf.jd > tr_ed)].time_min.to_numpy(),
                                      s_rf[(s_rf.jd < tr_st) | (s_rf.jd > tr_ed)][mg[idx]].to_numpy(),
                                      t=q_knots, s=1)

        yfit = interpolate.BSpline(tt, c, k)(s_rf.time_min.to_numpy())
        s_rf[cmg[idx]] = s_rf[mg[idx]] / yfit

        # plt.subplot(2, 1, 1)
        # plt.scatter(s_rf.time_min, s_rf[mg[idx]], c='k')
        # plt.plot(s_rf.time_min, yfit, c='b')

        # plt.subplot(2, 1, 2)
        # plt.plot(s_rf.time_min, s_rf[cmg[idx]], c=clr[(idx * 2) + 1],
        #         label='Reflection ' + Configuration.WAVELENGTHS_REFLECTION[idx])
        # plt.legend()
        # plt.show()

#for idx, row in s_tr.iterrows():
#    plt.scatter(Configuration.WAVELENGTHS_TRANSMISSION_NUMS, row[cmg] / s_tr[cmg].median(), c='k')
#    plt.scatter(Configuration.WAVELENGTHS_REFLECTION_NUMS, s_rf.loc[idx, cmg[:-1]] / s_rf[cmg[:-1]].median(), c='r')
#    plt.show()
plt.figure(figsize=(9, 6))
spect_660 = np.interp(660, spectra.wve_vc, spectra.titan)
for idx, m in enumerate(cmg):

    s_tr[clrmg[idx]] = s_tr[m] / s_tr['aper_c3']

    plt.scatter(s_tr['time_min'], (s_tr[clrmg[idx]] / s_tr[clrmg[idx]].median()) - idx / 2,
                c=clr[idx * 2], label=Configuration.WAVELENGTHS_TRANSMISSION[idx])
    plt.plot(s_tr['time_min'], (s_tr[clrmg[idx]] / s_tr[clrmg[idx]].median()) - idx / 2,
             c=clr[idx * 2], linewidth=2)

    if idx < 7:
        s_rf[clrmg[idx]] = s_rf[m] / s_tr['aper_c3']

        plt.scatter(s_rf['time_min'], (s_rf[clrmg[idx]] / s_rf[clrmg[idx]].median()) - (idx + 0.5) / 2,
                    c=clr[(idx * 2) + 1], label=Configuration.WAVELENGTHS_REFLECTION[idx])
        plt.plot(s_rf['time_min'], (s_rf[clrmg[idx]] / s_rf[clrmg[idx]].median()) - (idx + 0.5) / 2,
                 c=clr[(idx * 2) + 1], linewidth=2)

plt.ylabel('Color plus Offset [$\lambda$ / 660nm]', fontsize=15)
plt.xlabel('Time Since First Exposure [mins]', fontsize=15)
plt.legend()
plt.legend(fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
# plt.show()
plt.savefig('titan_photometry.png', bbox_inches='tight', pad_inches=0)
plt.close()

plt.figure(figsize=(9, 6))
spectrum = np.zeros(15)
for idx, m in enumerate(cmg):

    spectrum[idx * 2] = s_tr[m].min()
    if idx < 7:
        spectrum[(idx * 2) + 1] = s_rf[m].min()

plt.scatter(np.flip(wvs), spectrum / spectrum[4] * 100, c='k')
plt.plot(np.flip(wvs), spectrum / spectrum[4] * 100, c='k', linewidth=2, alpha=0.2)
# plt.scatter(wvs, spect / spect[10] * 100, c='r', label='Previous Measurement')
# plt.plot(spectra.wve_vc, spectra.titan / spect[10] * 100, c='r', linewidth=2, alpha=0.2)
# plt.scatter(wvs, spect2 / spect2[10] * 100, c='b', label='Previous Measurement')
# plt.plot(spectra2.wve/10, spectra2.flx / spect2[10] * 100, c='b', linewidth=2, alpha=0.2)
plt.scatter(wvs, spect3 / spect3[10] * 100, c='red', label='Neff+1984 Titan Flux')
plt.plot(spectra3.wve/10, spectra3.flx / spect3[10] * 100, c='red', linewidth=2, alpha=0.2)
plt.ylabel('Relative Occultation Depth [%, $\lambda$ / 660nm]', fontsize=15)
plt.xlabel('Wavelength [nm]', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
# plt.show()
plt.savefig('titan_spectrum.png', bbox_inches='tight', pad_inches=0)
plt.close()
