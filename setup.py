from distutils.core import setup, Extension

phred = Extension('phred', sources=['src/phred.c'], libraries=['m'])

setup(name='hoobari',
      packages=['hoobari'],
      package_dir= { 'hoobari':'src' },
      install_requires=['sqlite3','vcf','pickle','numpy','pandas','pysam'],
      ext_modules=[phred],
      version='1.0',
      author='Tom Rabinowitz',
      url='https://github.com/nshomron/hoobari',
      description='Bayesian fetal genotyping based on maternal cfDNA together with parental sequencing data')
