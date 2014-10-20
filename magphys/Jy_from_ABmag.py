from barak.sed import mag2Jy, get_bands

# AB MAGS

# they are galactic extinction corrected.

# Also IGM absorption corrected (u = -0.35, g = -0.08)

u =  26.25, 0.2
g =  26.02, 0.1
r =  25.9, 0.1
Ks = 25.1, 0.3

for val,sig in (u,g,r,Ks):
    jy = mag2Jy(val)
    jyhi = mag2Jy(val - sig)
    jylo = mag2Jy(val + sig)
    assert jyhi > 0
    assert jylo > 0
    jysig = 0.5 * (jyhi - jylo)
    print '  {:.5g}  {:.5g}'.format(jy, jysig),

print ''
print '24.5:', '{:.5g}'.format(mag2Jy(24.5))
