import distutils
from distutils.core import setup
import glob

bin_files = glob.glob("bin/*.py") 

# The main call
setup(name='wavg',
      version ='0.6',
      license = "GPL",
      description = "computes a WAVG[ERR]_MAG_PSF,  SPREAD[ERR]_MODEL quantities for each object in the coadd-catalog.",
      author = "Brian Yanny",
      author_email = "yanny@fnal.gov",
      packages = ['wavg'],
      package_dir = {'': 'python'},
      scripts = bin_files,
      data_files=[('ups',['ups/wavg.table'])],
      )

