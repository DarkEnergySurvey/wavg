#!/usr/bin/env python

import sys
import fitsio
import numpy
import subprocess
import math
import argparse

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
   

def readcoaddcat(filename,elems,hduforcoaddcat,urall,uraur,udecll,udecur):

  f=fitsio.FITS(filename,namemode='r')
  #This may have to be changed from HDU=1 to HDU=2
  t=f[int(hduforcoaddcat)].read()
 
  cd = {}
  for e in elems:
    cd[e] = t[e]


  cd['RA']=cd['ALPHAWIN_J2000']
  cd['NEPOCHS']=[0]*len(cd['RA'])
  cd['WAVG_MAG_PSF']=t['MAG_PSF']+0.0
  cd['WAVG_MAGERR_PSF']=t['MAGERR_PSF']+0.0
  cd['WAVG_MAGRMS_PSF']=t['MAGERR_PSF']+0.0
  cd['WAVG_SPREAD_MODEL']=t['SPREAD_MODEL']+0.0
  cd['WAVG_SPREADERR_MODEL']=t['SPREADERR_MODEL']+0.0
  cd['DUP']=[0]*len(cd['RA'])

  cindex=cd['ALPHAWIN_J2000'].argsort()

  #don't forget about the RA> 350 -360 wraparound

            
  coaddcrossra0 = 0
  if (float(urall) > -1 and float(urall) > 350 and float(uraur) < 10):
    coaddcrossra0 = 1
    uurall = float(urall)-360.0
    uuraur = float(uraur)
    uudecll = float(udecll)
    uudecur = float(udecur)
  else:
    uurall = float(urall)
    uuraur = float(uraur)
    uudecll = float(udecll)
    uudecur = float(udecur)

  if coaddcrossra0 == 1:
        for ii in range(0,len(cd['RA'])):
	  if cd['ALPHAWIN_J2000'][ii] >= 350:
		cd['RA'][ii] = cd['ALPHAWIN_J2000'][ii] - 360.0

  if (float(urall) > -1):
   lll=len(cd['RA'])
   for ee in range(0,lll):
     if (cd['RA'][ee] < uurall or cd['RA'][ee] > uuraur or cd['DELTAWIN_J2000'][ee] < uudecll or cd['DELTAWIN_J2000'][ee] > uudecur):
	cd['DUP'][ee] = 1


  cindex=cd['RA'].argsort()

  sortedra=cd['RA'][cindex]
  sorteddec=cd['DELTAWIN_J2000'][cindex]

  #print sortedra[0],sorteddec[0],sortedra[100],sorteddec[100],sortedra[1000],sorteddec[1000]
 
  #print "len nepochs:",len(cd['NEPOCHS']),len(cd['WAVG_SPREAD_MODEL'])

  return (cd,cindex)

def readlists(lists,sublistdirprefix,elems,hdrelems,urall,uraur,udecll,udecur):
	od={}
	for e in elems+hdrelems:
	 od[e] = []

	od['FILENAME'] = []
	od['MAG_ZERO'] = []
	od['COADD_NUMBER'] = []
	od['DUP'] = []

        alllists=lists.split(',')
        for list in alllists:
	 with open(list) as f:
	  for line in f:
	   #print line

	   if sublistdirprefix == 'none' or sublistdirprefix == 'None' or sublistdirprefix == 'NONE':
  	     getname=line.split(',')[0]
	     filename=getname.split('/')[-1]
	     zeropoint=float(line.split(',')[1])
	   else:
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
	     getname=sublistdirprefix+'/'+filename

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
	   od['DUP'] += [0]*sizetable

	   #print ra,dec,flux_psf
           #print ra[0],dec[0],flux_psf[0],expnum,ccdnum,ra.size,zeropoint,aband

	 #for e in elems+hdrelems:
	  #print e,len(od[e])
	 #print 'MAG_ZERO',len(od['MAG_ZERO'])
	 #print 'COADD_NUMBER',len(od['COADD_NUMBER'])

          nra = numpy.array(od['ALPHAWIN_J2000'])
	  oindex = nra.argsort()

	  crossra0 = 0
	  if (float(urall) > -1 and float(urall) > 350 and float(uraur) < 10):
            uurall = float(urall)-360.0
            uuraur = float(uraur)
            uudecll = float(udecll)
	    uudecur = float(udecur)
            crossra0 = 1
          else:
            uurall = float(urall)
            uuraur = float(uraur)
            uudecll = float(udecll)
	    uudecur = float(udecur)

	  #don't forget about the RA-360 for RA>350 wraparound.
  	  if crossra0 == 1:
             for ii in range(0,len(od['RA'])):
	       if od['ALPHAWIN_J2000'][ii] >= 350:
		 od['RA'][ii] = od['ALPHAWIN_J2000'][ii] - 360.0

	     #od['RA'] = od['ALPHAWIN_J2000'] if od['ALPHAWIN_J2000'] < 350 else od['ALPHAWIN_J2000'] - 360.0


	  if (float(urall) > -1):
	   lll=len(od['RA'])
           for ee in range(0,lll):
	     if (od['RA'][ee] < uurall or od['RA'][ee] > uuraur or od['DELTAWIN_J2000'][ee] < uudecll or od['DELTAWIN_J2000'][ee] > uudecur):
		od['DUP'][ee] = 1
             
	
          nra = numpy.array(od['RA'])

	  oindex = nra.argsort()

	  #print od['RA'][oindex[0]],od['RA'][oindex[1]],od['RA'][oindex[100000]]
	  #print onumber
    

	return (od,oindex)


def solvsort(cd,cindex,od,oindex,refmag0,matchradius):

  #takes two ra-sorted lists runs down the first list, flags matches to all things in the second within 1''

  totmatches=0
  sizec=len(cindex)
  sizeo=len(oindex)
  print "Size of Coadd Catalog:",sizec,"Total number of SE objects in lists:",sizeo
  topo=0
  curc=0
  tol=float(matchradius)/3600.0
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
    oneset['MAGRMS_PSF'] = []
    oneset['SPREAD_MODEL'] = []
    oneset['SPREADERR_MODEL'] = []
    oneset['MAG_ZERO'] = []

    while 1:

      if curo >= sizeo:
	break

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
	  oneset['MAGRMS_PSF'] += [od['MAGERR_PSF'][oindex[curo]]]
	  oneset['SPREAD_MODEL'] += [od['SPREAD_MODEL'][oindex[curo]]]
	  oneset['SPREADERR_MODEL'] += [od['SPREADERR_MODEL'][oindex[curo]]]
	  oneset['MAG_ZERO'] += [od['MAG_ZERO'][oindex[curo]]]

      curo += 1
 
    #here's where we compute the wavg numbers with this for this coadd object current set.
    (n,wavg_mag_psf,wavg_magerr_psf,wavg_magrms_psf,wavg_spread_model,wavg_spreaderr_model) =onebanderrs(oneset,refmag0)
    if n>0:
      #print '0',cd['MAG_PSF'][cindex[curc]],cd['MAGERR_PSF'][cindex[curc]], cd['SPREAD_MODEL'][cindex[curc]],cd['SPREADERR_MODEL'][cindex[curc]]
      cd['WAVG_MAG_PSF'][cindex[curc]] = wavg_mag_psf
      cd['WAVG_MAGERR_PSF'][cindex[curc]] = wavg_magerr_psf
      cd['WAVG_MAGRMS_PSF'][cindex[curc]] = wavg_magrms_psf
      cd['WAVG_SPREAD_MODEL'][cindex[curc]] = wavg_spread_model
      cd['WAVG_SPREADERR_MODEL'][cindex[curc]] = wavg_spreaderr_model
    else:
      cd['WAVG_MAG_PSF'][cindex[curc]] = -99
      cd['WAVG_MAGERR_PSF'][cindex[curc]] = -99
      cd['WAVG_MAGRMS_PSF'][cindex[curc]] = -99
      cd['WAVG_SPREAD_MODEL'][cindex[curc]] = -99
      cd['WAVG_SPREADERR_MODEL'][cindex[curc]] = -99
    curc += 1
      
  print "Number of SE objects matched to a Coadd object:",totmatches


def joinup(cd,cindex,od,oindex,outcoaddname,outocname,coaddcat):

  sizec=len(cindex)
  sizeo=len(oindex)
  #print sizec,sizeo

  ne=numpy.array(cd['NEPOCHS'])
  print "Distribution of NEPOCHS"
  print "NEPOCHs, Count(Coadd Objects with NEPOCHs)"

  for ii in range(0,21):
   print ii,ne[numpy.where(ne==ii)].size

  print ">20",ne[numpy.where(ne>10)].size

  #write out the table as a FITS file in order of NUMBER in the COADD catalog

  nrows=cd['NUMBER'].size
  band=od['BAND'][0]

  outdata=numpy.zeros(nrows,dtype=[('NUMBER','i4'),('NEPOCHS','i4'),('WAVG_MAG_PSF','f4'),('WAVG_MAGERR_PSF','f4'),('WAVG_MAGRMS_PSF','f4'),('WAVG_SPREAD_MODEL','f4'),('WAVG_SPREADERR_MODEL','f4'),('DUP','i4')])
  outdata['NUMBER']=cd['NUMBER']
  outdata['NEPOCHS']=cd['NEPOCHS']
  outdata['WAVG_MAG_PSF']=cd['WAVG_MAG_PSF']
  outdata['WAVG_MAGERR_PSF']=cd['WAVG_MAGERR_PSF']
  outdata['WAVG_MAGRMS_PSF']=cd['WAVG_MAGRMS_PSF']
  outdata['WAVG_SPREAD_MODEL']=cd['WAVG_SPREAD_MODEL']
  outdata['WAVG_SPREADERR_MODEL']=cd['WAVG_SPREADERR_MODEL']
  outdata['DUP']=cd['DUP']

  hlist=[{'name':'BAND','value':band}]
  fhdr=fitsio.FITSHDR(hlist)

  fitsio.write(outcoaddname,outdata,header=fhdr,extname='OBJECTS',clobber=True)

  #writeout the oc matching file:
  ne=numpy.array(od['COADD_NUMBER'])
  necid=ne[numpy.where(ne!=0)]
  no=numpy.array(od['NUMBER'])
  neoid=no[numpy.where(ne!=0)]
  nf=numpy.array(od['FILENAME'])
  nfid=nf[numpy.where(ne!=0)]
  ndup=numpy.array(od['DUP'])
  ndupid=ndup[numpy.where(ne!=0)]

  print "SE Objects with non-null coadd_number:",ne[numpy.where(ne!=0)].size

  nocrows=necid.size
  outocdata=numpy.zeros(nocrows,dtype=[('COADD_FILENAME','a72'),('NUMBER','i4'),('FILENAME','a60'),('SE_OBJECT_NUMBER','i4'),('DUP','i4')])
  outocdata['COADD_FILENAME']=[coaddcat]*nocrows
  outocdata['NUMBER']=necid
  outocdata['FILENAME']=nfid
  outocdata['SE_OBJECT_NUMBER']=neoid
  outocdata['DUP']=ndupid
  fitsio.write(outocname,outocdata,extname='OBJECTS',header=fhdr,clobber=True)

def onebanderrs(workdict,refmag0):

 m=[]
 e=[]
 w=[]
 sm=[]
 sme=[]
 smw=[]
 c=0.001	#error floor
 #c=0.0
 smc=0.000001
 #smc=0.0
 cnt=0
 weight_sum = 0.0
 sum_weights = 0.0
 sumsq_weights = 0.0

 smweight_sum = 0.0
 smsum_weights = 0.0
 smsumsq_weights = 0.0

 error_sum = 0.0
 error_smsum = 0.0

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

     error_sum += epsf*epsf

     sm.append(sm0)
     sme.append(sme0)
     smw.append(smw0)
     smweight_sum += smw0*sm0
     smsum_weights += smw0
     smsumsq_weights += smw0*smw0

     error_smsum += sme0*sme0

     cnt += 1
 
 n=cnt
 #print cnt,n
 if n==1:
  return (1,mpsf,epsf,epsf,sm0,sme0)
 if n>0:
  wavg = weight_sum/sum_weights
  weighted_devsq = 0
  weighted_sigma = numpy.sqrt(error_sum)/float(n)

  smwavg = smweight_sum/smsum_weights
  smweighted_devsq = 0
  smweighted_sigma = numpy.sqrt(error_smsum)/float(n)

  for i in range(0,n):
   weighted_devsq += (m[i]-wavg)*(m[i]-wavg)*w[i]
   smweighted_devsq += (sm[i]-smwavg)*(sm[i]-smwavg)*smw[i]

  weighted_rms = math.sqrt(weighted_devsq/(sum_weights-(sumsq_weights/sum_weights)))
  smweighted_rms = math.sqrt(smweighted_devsq/(smsum_weights-(smsumsq_weights/smsum_weights)))
 
  #print "n,weighted average,weighted_sigma:",n,wavg,weighted_sigma

  return  (n,wavg,weighted_sigma,weighted_rms,smwavg,smweighted_sigma)

 else:
  return  (0,0,0,0,0,0)

#END OF ONEBANDERRS

parser=argparse.ArgumentParser()
parser.add_argument("incoaddcat",help="Filename of the input coadd fits catalog (with no path). Also used as part of unique coadd object identifier.")
parser.add_argument("catlists",help="List(s) of single epoch catalogs overlapping the tile.")
parser.add_argument("outwavgcat",help="filename of the output fits file with wavg quantities.")
parser.add_argument("outoclinkcat",help="filename of the output fits file with coadd_object ids for each object id that matches a coadd_object.")
parser.add_argument("--urall",help="unique ra lower left (deg)",default=-1.0)
parser.add_argument("--uraur",help="unique ra upper right (deg)",default=-1.0)
parser.add_argument("--udecll",help="unique dec lower left (deg)",default=-1.0)
parser.add_argument("--udecur",help="unique dec upper right (deg)",default=-1.0)
parser.add_argument("--dircoaddcat",help="directory containing the coadd catalog file (defaults to current dir .)",default=".")
parser.add_argument("--sublistdirprefix",help="perform a substitution on the list directory changing red/ to sublistdirprefix/cat/ and immasked --> red-fullcat containing the coadd catalog file (defaults to no substitution in list)",default="none")
parser.add_argument("--refmag0",help="sextractor single epoch reference magnitude (default 25.0)",default=25.0)
parser.add_argument("--matchradius",help="match radius in arcsec for single obj, coadd obj matching (default 0.7)", default=0.7)
parser.add_argument("--hduforcoaddcat",help="hdu for coadd catalog (default 1)",default=1)
args=parser.parse_args()

#fetch of catalogs not used now -- assumed to be present in working area
#grabredfullcats(listname)

elems = ['NUMBER','ALPHAWIN_J2000','DELTAWIN_J2000','SPREAD_MODEL','SPREADERR_MODEL','FLAGS','IMAFLAGS_ISO','MAG_PSF','MAGERR_PSF']

hdrelems = ['BAND','CCDNUM','EXPNUM']

(cd,cindex) = readcoaddcat(args.dircoaddcat+'/'+args.incoaddcat,elems,args.hduforcoaddcat,args.urall,args.uraur,args.udecll,args.udecur)

(od,oindex) = readlists(args.catlists,args.sublistdirprefix,elems,hdrelems,args.urall,args.uraur,args.udecll,args.udecur)

solvsort(cd,cindex,od,oindex,args.refmag0,args.matchradius)
joinup(cd,cindex,od,oindex,args.outwavgcat,args.outoclinkcat,args.incoaddcat)

#done
