#gain=0.45
# r band
seeing=0.61
detect=1.5
analysis=1.5
deblend=0.00001
btype=AUTO    #MANUAL
bsize=128
# aperture diameter n pixels
aperturepix=11
# 6 arcsec

# a catalogue for each band for photometry purposes, then catalogues
# detect in R and measureing photmetry in u and g.

############# 

ref=r

band=g

prefix1=$ref
prefix2=$band
image=$prefix1.fits,$prefix2.fits
catalogue=$band$ref.cat
sex $image -c conf/phot.sex -CATALOG_NAME $catalogue -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME conf/phot.param -DETECT_THRESH $detect -ANALYSIS_THRESH $analysis -FILTER_NAME conf/default.conv -STARNNW_NAME conf/default.nnw -VERBOSE_TYPE NORMAL -DEBLEND_MINCONT $deblend -BACK_TYPE $btype -BACK_SIZE $bsize -MEMORY_BUFSIZE 4096 -SEEING_FWHM $seeing -PHOT_APERTURES $aperturepix
#-WEIGHT_TYPE NONE 

mv check.fits check_$prefix1$prefix2.fits

band=u

prefix1=$ref
prefix2=$band
image=$prefix1.fits,$prefix2.fits
catalogue=$band$ref.cat
sex $image -c conf/phot.sex -CATALOG_NAME $catalogue -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME conf/phot.param -DETECT_THRESH $detect -ANALYSIS_THRESH $analysis -FILTER_NAME conf/default.conv -STARNNW_NAME conf/default.nnw -VERBOSE_TYPE NORMAL -DEBLEND_MINCONT $deblend -BACK_TYPE $btype -BACK_SIZE $bsize -MEMORY_BUFSIZE 4096 -SEEING_FWHM $seeing -PHOT_APERTURES $aperturepix
#-WEIGHT_TYPE NONE 

mv check.fits check_$prefix1$prefix2.fits

# single catalogues

################# 
# band=u

# prefix=$band
# image=$prefix.fits
# catalogue=$band.cat
# sex $image -c conf/phot.sex -CATALOG_NAME $catalogue -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME conf/phot.param -DETECT_THRESH $detect -ANALYSIS_THRESH $analysis -FILTER_NAME conf/default.conv -STARNNW_NAME conf/default.nnw -VERBOSE_TYPE NORMAL -DEBLEND_MINCONT $deblend -BACK_TYPE AUTO -MEMORY_BUFSIZE 4096 -SEEING_FWHM $seeing
# #-PHOT_APERTURES $aperturepix
# #-WEIGHT_TYPE NONE 


# band=g

# prefix=$band
# image=$prefix.fits
# catalogue=$band.cat
# sex $image -c conf/phot.sex -CATALOG_NAME $catalogue -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME conf/phot.param -DETECT_THRESH $detect -ANALYSIS_THRESH $analysis -FILTER_NAME conf/default.conv -STARNNW_NAME conf/default.nnw -VERBOSE_TYPE NORMAL -DEBLEND_MINCONT $deblend -BACK_TYPE AUTO -MEMORY_BUFSIZE 4096 -SEEING_FWHM $seeing
# #-PHOT_APERTURES $aperturepix
# #-WEIGHT_TYPE NONE 


band=r

prefix=$band
image=$prefix.fits
catalogue=$band.cat
sex $image -c conf/phot.sex -CATALOG_NAME $catalogue -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME conf/phot.param -DETECT_THRESH $detect -ANALYSIS_THRESH $analysis -FILTER_NAME conf/default.conv -STARNNW_NAME conf/default.nnw -VERBOSE_TYPE NORMAL -DEBLEND_MINCONT $deblend -BACK_TYPE $btype -BACK_SIZE $bsize -MEMORY_BUFSIZE 4096 -SEEING_FWHM $seeing -PHOT_APERTURES $aperturepix
#-WEIGHT_TYPE NONE 

mv check.fits check_$prefix.fits

sex2DS9reg $catalogue $band.reg
