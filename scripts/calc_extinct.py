# The use of IGMtransmission should be acknowledged using a statement
# like This paper computed IGM transmission values using
# IGMtransmissiona (Harrison, Meiksin & Stock 2011), based on the
# transmission curves of Meiksin (2006).  aAvailable for download from
# http://code.google.com/p/igmtransmission If the LLS distribution
# from Inoue & Iwata (2008) is used, a reference to their work should
# be added as well.

from barak.io import readtxt
import numpy as np

# z=2.76, Inoue LLS code. Diffuse IGM normalisation 0.07553

T = readtxt('averageTransmission.dat', names='wa,tr')
WA = np.arange(3180, 8000, 1.)
TR = np.interp(WA, T.wa, T.tr)
TR[WA > 4543.] = 1

from barak.sed import get_bands
u,g = get_bands('sdss', ['u','g'])

utr = np.interp(WA, u.wa, u.tr)
av_utr = (utr * TR).sum() / utr.sum()
print 'u bandpass-weighted transmission', av_utr
print 'mag extinct', -2.5*np.log10(av_utr)

gtr = np.interp(WA, g.wa, g.tr)
av_gtr = (gtr * TR).sum() / gtr.sum()
print 'g bandpass-weighted transmission', av_gtr
print 'mag extinct', -2.5*np.log10(av_gtr)
