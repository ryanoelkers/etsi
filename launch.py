import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression


def custom_round(x, base=5):
    return int(base * round(float(x)/base))


date = '2022-06-20'
s20 = pd.read_csv("C:\\Users\\barristan\\Development\\etsi\\analysis\\Nicole\\0620\\Nicole_transmission_00.lc",
                  sep=' ', header=0)
c20 = pd.read_csv("C:\\Users\\barristan\\Development\\etsi\\analysis\\Nicole\\0620\\Nicole_transmission_01.lc",
                  sep=' ', header=0)

date = '2022-06-22'
s22 = pd.read_csv("C:\\Users\\barristan\\Development\\etsi\\analysis\\Nicole\\0622\\Nicole_transmission_00.lc",
                  sep=' ', header=0)
c22 = pd.read_csv("C:\\Users\\barristan\\Development\\etsi\\analysis\\Nicole\\0622\\Nicole_transmission_01.lc",
                  sep=' ', header=0)

s22['time_min'] = s22.apply(lambda x: (x.jd - s22.jd[0]) * 24. * 60., axis=1)
c22['time_min'] = c22.apply(lambda x: (x.jd - c22.jd[0]) * 24. * 60., axis=1)

s20['time_min'] = s20.apply(lambda x: (x.jd - s20.jd[0]) * 24. * 60., axis=1)
c20['time_min'] = c20.apply(lambda x: (x.jd - c20.jd[0]) * 24. * 60., axis=1)

s22 = s22[s22.time_min < 350].reset_index(drop=True)
c22 = c22[c22.time_min < 350].reset_index(drop=True)

s22 = s22.interpolate()
c22 = c22.interpolate()

s17 = s17[(c17.time_min > 5) & (s17.time_min < 350)].reset_index(drop=True)
c17 = c17[(c17.time_min > 5) & (c17.time_min < 350)].reset_index(drop=True)
s17 = s17.interpolate()
c17 = c17.interpolate()

mg = ['mag_p1', 'mag_p2', 'mag_p3', 'mag_p4', 'mag_p5', 'mag_p6', 'mag_p7', 'mag_p8']
pos_x = ['x_p1', 'x_p2', 'x_p3', 'x_p4', 'x_p5', 'x_p6', 'x_p7', 'x_p8']
pos_y = ['y_p1', 'y_p2', 'y_p3', 'y_p4', 'y_p5', 'y_p6', 'y_p7', 'y_p8']

s22['binned'] = s22.apply(lambda x: custom_round(x.time_min), axis=1)
c22['binned'] = c22.apply(lambda x: custom_round(x.time_min), axis=1)
s17['binned'] = s17.apply(lambda x: custom_round(x.time_min), axis=1)
c17['binned'] = c17.apply(lambda x: custom_round(x.time_min), axis=1)
clr = ['maroon', 'red', 'salmon', 'darkorange', 'gold', 'darkgreen', 'blue', 'purple']

for idx, m in enumerate(mg):

    # do the June 17 night
    z = c17[m].rolling(300, min_periods=1, center=True).median().to_numpy().reshape(-1, 1)
    y = s17[m].to_numpy()
    x = c17[m].to_numpy().reshape(-1, 1)

    reg = LinearRegression().fit(z, y)
    yhat17 = reg.predict(x)

    s17[m] = s17[m]-(yhat17 - np.median(yhat17) + np.median(s17[m]))
    s_bin_17 = s17.groupby('binned').agg({m: 'median'}).reset_index()

    # do the June 22 night
    z = c22[m].rolling(300, min_periods=1, center=True).median().to_numpy().reshape(-1, 1)
    y = s22[m].to_numpy()
    x = c22[m].to_numpy().reshape(-1, 1)

    reg = LinearRegression().fit(z, y)
    yhat22 = reg.predict(x)

    s22[m] = s22[m]-(yhat22 - np.median(yhat22) + np.median(s22[m]))
    s_bin_22 = s22.groupby('binned').agg({m: 'median'}).reset_index()

    plt.scatter(s_bin_17.binned, s_bin_17[m], label='Global', marker='x', c=clr[idx])
    plt.scatter(s_bin_22.binned, s_bin_22[m], label='Global+Local', marker='o', c=clr[idx])
    plt.xlabel('Time from Observation Start [mins]')
    plt.ylabel('Relative Magnitude')

plt.legend()
plt.ylim([0.04, -0.04])
plt.show()


s_bin_22 = s22.groupby('binned').agg({'mag_p3': 'median', 'mag_p2': 'median', 'mag_p7': 'median'}).reset_index()
s_bin_17 = s17.groupby('binned').agg({'mag_p3': 'median', 'mag_p2': 'median', 'mag_p7': 'median'}).reset_index()
res1 = np.polyfit(s_bin_17.binned, s_bin_17['mag_p3'] - s_bin_17['mag_p2'], 1)
plt.scatter(s_bin_17.binned, s_bin_17['mag_p3'] - s_bin_17['mag_p2'], label='Global', marker='x', c='blue')
plt.plot(s_bin_17.binned, res1[0] * s_bin_17.binned + res1[1], c='b')
bb = s_bin_22['mag_p3'] - s_bin_22['mag_p2']
res2 = np.polyfit(s_bin_22[(bb > -0.04) & (bb < 0.04)].binned,
                  s_bin_22[(bb > -0.04) & (bb < 0.04)]['mag_p3'] - s_bin_22[(bb > -0.04) & (bb < 0.04)]['mag_p2'], 1)
plt.scatter(s_bin_22.binned, s_bin_22['mag_p3'] - s_bin_22['mag_p2'], label='Global+Local', marker='o', c='red')
plt.plot(s_bin_22.binned, res2[0] * s_bin_22.binned + res2[1], c='r')
plt.xlabel('Time from Observation Start [mins]')
plt.ylabel('[435nm - 763nm] [mag]')

plt.legend()
plt.ylim([0.04, -0.04])
plt.show()

comp = s22.groupby('binned').agg({'mag_p2': 'median'}).reset_index()
for idx, m in enumerate(mg):
    s_bin = s22.groupby('binned').agg({m: 'median'}).reset_index()
    plt.scatter(s_bin.binned, s_bin[m] - comp['mag_p2'], label=m)
    plt.legend()
    plt.gca().invert_yaxis()
    plt.show()
    print('hold')


s17 = pd.read_csv("C:\\Users\\barristan\\Desktop\\XO-1_transmission_01_17.lc", sep=' ', header=0)
c17 = pd.read_csv("C:\\Users\\barristan\\Desktop\\XO-1_transmission_00_17.lc", sep=' ', header=0)

s17['time_min'] = s17.apply(lambda x: (x.jd - s17.jd[0]) * 24. * 60., axis=1)
c17['time_min'] = c17.apply(lambda x: (x.jd - c17.jd[0]) * 24. * 60., axis=1)

