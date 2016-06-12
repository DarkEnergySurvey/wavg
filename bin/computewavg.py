#!/usr/bin/env python

import sys
import fitsio
import numpy
import subprocess
import math

def Usage():
	sys.exit("computewavg list tile coaddcat outcoaddname outocname refmagzero coaddcatinname")

def grabredfullcats(list):
 prefix='https://desar2.cosmology.illinois.edu/DESFiles/Prodalpha/archive/ACT/finalcut/Y2T5-2112'
 with open(list) as f:
  for line in f:
   #print line
   a=line.split(',')[0]
   b=a.split('/')[1]
   c=b.split('_')
   expnumstring=c[0]
   band=c[1]
   ccdstring=c[2]
   rpstring=c[3]
   pstring=c[3].split('p')[1]
   getname=prefix+'/'+expnumstring+'/p'+pstring+'/cat/'+expnumstring+'_'+band+'_'+ccdstring+'_'+rpstring+'_red-fullcat.fits'
   subprocess.check_output(['wget',getname])
   

def readcoaddcat(filename,elems):

  f=fitsio.FITS(filename,namemode='r')
  #This may have to be changed from HDU=2 to HDU=1
  t=f[2].read()
 
  cd = {}
  for e in elems:
    cd[e] = t[e]

  cd['RA']=cd['ALPHAWIN_J2000']
  cd['NEPOCHS']=[0]*len(cd['RA'])
  cd['WAVG_MAG_PSF']=t['MAG_PSF']
  cd['WAVG_MAGERR_PSF']=t['MAGERR_PSF']
  cd['WAVG_SPREAD_MODEL']=t['SPREAD_MODEL']
  cd['WAVG_SPREADERR_MODEL']=t['SPREADERR_MODEL']

  cindex=cd['ALPHAWIN_J2000'].argsort()

  #don't forget about the RA> 350 -360 wraparound

  if cd['ALPHAWIN_J2000'][cindex[0]] < 10 and cd['ALPHAWIN_J2000'][cindex[len(cindex)-1]] > 350:
	cd['RA'] = cd['ALPHAWIN_J2000'] if cd['ALPHAWIN_J2000'] < 350 else cd['ALPHAWIN_J2000'] - 360.0

  cindex=cd['RA'].argsort()

  sortedra=cd['RA'][cindex]
  sorteddec=cd['DELTAWIN_J2000'][cindex]

  #print sortedra[0],sorteddec[0],sortedra[100],sorteddec[100],sortedra[1000],sorteddec[1000]
 
  #print "len nepochs:",len(cd['NEPOCHS']),len(cd['WAVG_SPREAD_MODEL'])

  return (cd,cindex)

def readlists(list,tile,elems,hdrelems):
	od={}
	for e in elems+hdrelems:
	 od[e] = []

	od['FILENAME'] = []
	od['MAG_ZERO'] = []
	od['COADD_NUMBER'] = []

	with open(list) as f:
	 for line in f:
	  #print line

  	  a=line.split('/')[1]
  	  #print a
  	  b=a.split(',')[0]
	  zeropoint=float(a.split(',')[1])
  	  #print b
          c=b.split('_')
          expnumstring=c[0]
          nband=c[1]
          ccdstring=c[2]
          rpstring=c[3]
          pstring=c[3].split('p')[1]
	  filename= expnumstring+'_'+nband+'_'+ccdstring+'_'+rpstring+'_red-fullcat.fits'
	  getname=tile+'/cat/'+filename

	  fits=fitsio.FITS(getname.strip(),namemode='r')

	  hdu=0
	  h=fits[hdu].read_header()

	  table = fits[hdu+2].read()

	  #these are arrays

	  sizetable = table['NUMBER'].size

	  for e in elems:
	   od[e] += table[e].tolist()

	  od['RA'] = od['ALPHAWIN_J2000']
 
	  for he in hdrelems:
	   od[he] += [h[he]]*sizetable
	  

	  od['FILENAME'] += [filename]*sizetable
	  od['MAG_ZERO'] += [zeropoint]*sizetable
	  od['COADD_NUMBER'] += [0]*sizetable

	  #print ra,dec,flux_psf
          #print ra[0],dec[0],flux_psf[0],expnum,ccdnum,ra.size,zeropoint,aband

	 #for e in elems+hdrelems:
	  #print e,len(od[e])
	 #print 'MAG_ZERO',len(od['MAG_ZERO'])
	 #print 'COADD_NUMBER',len(od['COADD_NUMBER'])

         nra = numpy.array(od['ALPHAWIN_J2000'])
	 oindex = nra.argsort()

	 #don't forget about the RA-360 for RA>350 wraparound.
  	 if od['ALPHAWIN_J2000'][oindex[0]] < 10 and od['ALPHAWIN_J2000'][oindex[len(oindex)-1]] > 350:
	    od['RA'] = od['ALPHAWIN_J2000'] if od['ALPHAWIN_J2000'] < 350 else od['ALPHAWIN_J2000'] - 360.0

         nra = numpy.array(od['RA'])
	 oindex = nra.argsort()

	 #print od['RA'][oindex[0]],od['RA'][oindex[1]],od['RA'][oindex[100000]]
	 #print onumber
    

	return (od,oindex)


def solvsort(cd,cindex,od,oindex,refmag0):

  #takes two ra-sorted lists runs down the first list, flags matches to all things in the second within 1''

  totmatches=0
  sizec=len(cindex)
  sizeo=len(oindex)
  print sizec,sizeo
  topo=0
  curc=0
  tol=0.7/3600.0
  tolsq=tol*tol
  tol2=3*tol
  dtr=3.14159/180

  
  while curc < sizec:
  #while curc < 10:

    cra=cd['RA'][cindex[curc]]
    cdec=cd['DELTAWIN_J2000'][cindex[curc]]
    cosd=math.cos(cdec*dtr)
    cnumber=cd['NUMBER'][cindex[curc]]

    #print curc,cra,cdec,cnumber,cd['FLUX_PSF'][cindex[curc]],cd['FLUXERR_PSF'][cindex[curc]]

    if topo >= sizeo:
	break

    curo=topo

    #initialize current set dict for one coadd object's matches:
    oneset={}
    oneset['FLAGS'] = []
    oneset['MAG_PSF'] = []
    oneset['MAGERR_PSF'] = []
    oneset['SPREAD_MODEL'] = []
    oneset['SPREADERR_MODEL'] = []
    oneset['MAG_ZERO'] = []

    while 1:

      ora=od['RA'][oindex[curo]]
      odec=od['DELTAWIN_J2000'][oindex[curo]]

      delra=cra-ora

      if delra > tol2:
        topo += 1
        curo += 1

      if delra < -tol2:
        break

      deltasq=(cra-ora)*(cra-ora)*cosd*cosd+(cdec-odec)*(cdec-odec)

      if deltasq < tolsq:
	if od['FLAGS'][oindex[curo]] < 4 and od['IMAFLAGS_ISO'][oindex[curo]] == 0 and od['MAG_PSF'][oindex[curo]] < 99 and od['COADD_NUMBER'][oindex[curo]] == 0 :
	  #got a match! clean flags and imaflags and such
	  #print "match:",cnumber,curc,curo,ora,odec,od['CCDNUM'][oindex[curo]],od['EXPNUM'][oindex[curo]],od['FLUX_PSF'][oindex[curo]],od['MAG_ZERO'][oindex[curo]]
	  totmatches += 1
          od['COADD_NUMBER'][oindex[curo]]=cnumber
	  cd['NEPOCHS'][cindex[curc]] += 1
	  oneset['FLAGS'] += [od['FLAGS'][oindex[curo]]]
	  oneset['MAG_PSF'] += [od['MAG_PSF'][oindex[curo]]]
	  oneset['MAGERR_PSF'] += [od['MAGERR_PSF'][oindex[curo]]]
	  oneset['SPREAD_MODEL'] += [od['SPREAD_MODEL'][oindex[curo]]]
	  oneset['SPREADERR_MODEL'] += [od['SPREADERR_MODEL'][oindex[curo]]]
	  oneset['MAG_ZERO'] += [od['MAG_ZERO'][oindex[curo]]]

      curo += 1
 
    #here's where we compute the wavg numbers with this for this coadd object current set.
    (n,wavg_mag_psf,wavg_magerr_psf,wavg_spread_model,wavg_spreaderr_model) =onebanderrs(oneset,refmag0)
    #print n,wavg_mag_psf,wavg_magerr_psf,wavg_spread_model,wavg_spreaderr_model
    if n>0:
      #print '0',cd['MAG_PSF'][cindex[curc]],cd['MAGERR_PSF'][cindex[curc]], cd['SPREAD_MODEL'][cindex[curc]],cd['SPREADERR_MODEL'][cindex[curc]]
      cd['WAVG_MAG_PSF'][cindex[curc]] = wavg_mag_psf
      cd['WAVG_MAGERR_PSF'][cindex[curc]] = wavg_magerr_psf
      cd['WAVG_SPREAD_MODEL'][cindex[curc]] = wavg_spread_model
      cd['WAVG_SPREADERR_MODEL'][cindex[curc]] = wavg_spreaderr_model
    curc += 1
      
  print "totmatches:",totmatches


def joinup(cd,cindex,od,oindex,outcoaddname,outocname,coaddcat):

  sizec=len(cindex)
  sizeo=len(oindex)
  print sizec,sizeo

  ne=numpy.array(cd['NEPOCHS'])
  for ii in range(0,11):
   print ii,ne[numpy.where(ne==ii)].size

  print ">10",ne[numpy.where(ne>10)].size


  #write out the table as a FITS file in order of NUMBER in the COADD catalog

  nrows=cd['NUMBER'].size

  outdata=numpy.zeros(nrows,dtype=[('NUMBER','i4'),('NEPOCHS','i4'),('WAVG_MAG_PSF','f4'),('WAVG_MAGERR_PSF','f4'),('WAVG_SPREAD_MODEL','f4'),('WAVG_SPREADERR_MODEL','f4')])
  outdata['NUMBER']=cd['NUMBER']
  outdata['NEPOCHS']=cd['NEPOCHS']
  outdata['WAVG_MAG_PSF']=cd['WAVG_MAG_PSF']
  outdata['WAVG_MAGERR_PSF']=cd['WAVG_MAGERR_PSF']
  outdata['WAVG_SPREAD_MODEL']=cd['WAVG_SPREAD_MODEL']
  outdata['WAVG_SPREADERR_MODEL']=cd['WAVG_SPREADERR_MODEL']
  fitsio.write(outcoaddname,outdata,clobber=True)

  #writeout the oc matching file:
  ne=numpy.array(od['COADD_NUMBER'])
  necid=ne[numpy.where(ne!=0)]
  no=numpy.array(od['NUMBER'])
  neoid=no[numpy.where(ne!=0)]
  nf=numpy.array(od['FILENAME'])
  nfid=nf[numpy.where(ne!=0)]

  print "non-zero coadd_number:",ne[numpy.where(ne!=0)].size

  nocrows=necid.size
  outocdata=numpy.zeros(nocrows,dtype=[('COADD_FILENAME','a72'),('COADD_NUMBER','i4'),('FILENAME','a60'),('NUMBER','i4')])
  outocdata['COADD_FILENAME']=[coaddcat]*nocrows
  outocdata['COADD_NUMBER']=necid
  outocdata['FILENAME']=nfid
  outocdata['NUMBER']=neoid
  fitsio.write(outocname,outocdata,clobber=True)

def onebanderrs(workdict,refmag0):

 m=[]
 e=[]
 w=[]
 sm=[]
 sme=[]
 smw=[]
 c=0.01	#error floor
 smc=0.00001
 cnt=0
 weight_sum = 0.0
 sum_weights = 0.0
 sumsq_weights = 0.0

 smweight_sum = 0.0
 smsum_weights = 0.0
 smsumsq_weights = 0.0

 #print 'here:',len(workdict['FLAGS'])
 flagarr=numpy.array(workdict['FLAGS'])
 mpsffarr=numpy.array(workdict['MAG_PSF'])
 zparr=numpy.array(workdict['MAG_ZERO'])
 mpsferrarr=numpy.array(workdict['MAGERR_PSF'])
 smarr=numpy.array(workdict['SPREAD_MODEL'])
 smerrarr=numpy.array(workdict['SPREADERR_MODEL'])
 for i in range (0,len(flagarr)):
    flags=int(flagarr[i])
    mpsf=float(mpsffarr[i])+float(zparr[i])-float(refmag0)
    epsf=float(mpsferrarr[i])
    wpsf=1.0/(epsf*epsf+c*c)
    sm0=float(smarr[i])
    sme0=float(smerrarr[i])
    #print 'i,mpsf,sm',i,mpsf,sm0
    smw0=1.0/(sme0*sme0+smc*smc)

    #print i,mpsf,epsf,sm0,sme0
    if mpsf < 99 and flags < 4:
     m.append(mpsf)
     e.append(epsf)
     w.append(wpsf)
     weight_sum += wpsf*mpsf
     sum_weights += wpsf
     sumsq_weights += wpsf*wpsf

     sm.append(sm0)
     sme.append(sme0)
     smw.append(smw0)
     smweight_sum += smw0*sm0
     smsum_weights += smw0
     smsumsq_weights += smw0*smw0

     cnt += 1
 
 n=cnt
 #print cnt,n
 if n==1:
  return (1,mpsf,epsf,sm0,sme0)
 if n>0:
  wavg = weight_sum/sum_weights
  weighted_devsq = 0

  smwavg = smweight_sum/smsum_weights
  smweighted_devsq = 0

  for i in range(0,n):
   weighted_devsq += (m[i]-wavg)*(m[i]-wavg)*w[i]
   smweighted_devsq += (sm[i]-smwavg)*(sm[i]-smwavg)*smw[i]

  weighted_sigma = math.sqrt(weighted_devsq/(sum_weights-(sumsq_weights/sum_weights)))
  smweighted_sigma = math.sqrt(smweighted_devsq/(smsum_weights-(smsumsq_weights/smsum_weights)))
 
  #print "n,weighted average,weighted_sigma:",n,wavg,weighted_sigma

  return  (n,wavg,weighted_sigma,smwavg,smweighted_sigma)

 else:
  return  (0,0,0,0,0)

#END OF ONEBANDERRS

def computewavg(listname,tile,coaddcat,outcoaddname,outocname,refmag0,coaddcatname):

  #grabredfullcats(listname)

  elems = ['NUMBER','ALPHAWIN_J2000','DELTAWIN_J2000','FLUX_PSF','FLUXERR_PSF','SPREAD_MODEL','SPREADERR_MODEL','FLAGS','IMAFLAGS_ISO','MAG_PSF','MAGERR_PSF']
  hdrelems = ['BAND','CCDNUM','EXPNUM']

  (cd,cindex) = readcoaddcat(coaddcat,elems)

  (od,oindex) = readlists(listname,tile,elems,hdrelems)

  solvsort(cd,cindex,od,oindex,refmag0)
  joinup(cd,cindex,od,oindex,outcoaddname,outocname,coaddcatname)

  #writeout()


if __name__ == "__main__":
  if len(sys.argv[1:]) < 7 or (sys.argv[1] == '-h'):
     Usage()
  computewavg(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])

