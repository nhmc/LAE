from scipy.optimize import curve_fit
from astropy.io import fits
import numpy as np
import pylab as plt

def Gauss2d(X, x0, y0, sigx, sigy):
    x,y = X
    return np.exp(-((x-x0)/sigx)**2 - ((y - y0)/sigy)**2)

D = fits.getdata('r.fits')

if 1:
    dx = dy = 200
    x0 = 520
    y0 = 350
if 0:
    dx = dy = 45
    x0 = 590
    y0 = 420


x1 = x0 + dx
y1 = y0 + dy

d = D[y0:y1,x0:x1]
plt.figure(1)
plt.clf()
ax = plt.gca()
ax.matshow(np.arcsinh(d),cmap=plt.cm.hot)

xc = 91.4 - (x0-520)
yc = 91.4 - (y0-350)

ax.plot(xc,yc, 'x', ms=10)

step = 60

X = np.meshgrid(np.arange(dy), np.arange(dx))

ax.axis((xc-step,xc+step,yc-step,yc+step))

def func(X, h1,h2,h3,h4,xc,yc,s1,s2,s3,s4):
    g = h1*Gauss2d(X,xc,yc,s1,s1) + \
           h2*Gauss2d(X,xc,yc,s2,s2) + \
           h3*Gauss2d(X,xc,yc,s3,s3) + \
           h4*Gauss2d(X,xc,yc,s4,s4)
    return g.ravel()

#def func(X, h1,xc,yc,s1,):
#    g = h1*Gauss2d(X,xc,yc,s1)
#    return g.ravel()

    
if 1:
    # h1,h2,h3,h4,xc,yc,s1,s2,s3,s4
    #popt =  2e2,1.8e3,1.1e4,1.6e5, 21.31, 21.28, 15, 7, 4.5, 2.5
    # popt = np.array([  1.36401522e+02,   1.38327111e+03,   1.61291810e+04,
    #      1.49072098e+05,   2.13102440e+01 - (x0 - 590),
    #                    2.12792483e+01 - (y0 - 420),
    #      1.68451571e+01,   7.96521960e+00,   4.35016052e+00,
    #      2.61810186e+00])
    popt = np.array([  1.36401522e+02,   1.28327111e+03,   1.55e4,
         1.49e5,   2.13102440e+01 - (x0 - 590),
                       2.12792483e+01 - (y0 - 420),
         1.68451571e+01,   7.96521960e+00,   4.35016052e+00,
         2.61810186e+00])
    plt.figure(2)
    plt.clf()
    ax = plt.gca()
    g = func(X, *popt).reshape(dy,dx)
    ax.matshow(np.arcsinh((d-g))[::-1, :],cmap=plt.cm.hot)
#ax.matshow(np.arcsinh(d-(g1+g2)),cmap=plt.cm.hot)


if 0:
    h1,h2,h3,h4,xc,yc,s1,s2,s3,s4 = (
        2e2,1.8e3,1.1e4,1.6e5, 21.31, 21.28,15, 7, 4.5, 2.5)
    p0 = (h1,h2,h3,h4,xc,yc,s1,s2,s3,s4)
    popt, pcov = curve_fit(func, X, d.ravel(), p0=p0)
    plt.figure(3)
    plt.clf()
    ax = plt.gca()
    g = func(X, *popt).reshape(dy,dx)
    ax.matshow(np.arcsinh((d-g))[::-1, :],cmap=plt.cm.hot)

plt.show()

fh = fits.open('r.fits')
D[y0:y1,x0:x1] = d-g
fh[0].data = D
fh.writeto('rsub1.fits', clobber=1)

