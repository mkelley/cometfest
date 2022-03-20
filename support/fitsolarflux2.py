from scipy import signal
from numpy import *
from numpy.polynomial.chebyshev import chebfit, chebval
from matplotlib.pyplot import *
import mskpy

# fit a polynomial to the solar spectrum at 1 AU

w, s = mskpy.calib.e490()  # W/m2/um
s *= w**2 / mskpy.constants.c * 1e26  # Jy
kern = mskpy.Gaussian(linspace(0, 1, 100), 0.5, 0.1)
kern /= kern.sum()
ss = signal.convolve(s, kern, mode='same')
# the smoothed spectrum is good until ~10 um
S = r_[ss[w < 9], s[w >= 9]]

wr = [0, 0.25, 0.3, 0.4, 0.6, 1.0, 3, 15, 1000]
C, T = zeros((2, len(wr) - 1))
Sfit = zeros(len(S))
for i in range(len(wr) - 1):
    j = (w > wr[i]) * (w <= wr[i + 1])
    T[i], C[i] = mskpy.planckfit(w[j], S[j], S[j] * 0.01, (5800.0, 1.0))[0]

    Sfit[j] = mskpy.Planck(w[j], T[i]) * C[i]

#fit = mskpy.planckfit(w, S, S * 0.01, (5770, 1.0))
#Sfit = mskpy.Planck(w, fit[0][0]) * fit[0][1]

figure(1)
clf()
plot(w, Sfit / S)
axhline(1, color='r', alpha=0.5)
setp(gca(), xscale='log', yscale='log', ylim=[0.1, 10])
draw()


