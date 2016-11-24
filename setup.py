
__authors__ = ["FM Surace, C Richter, D Novikov, M Sanchez del Rio"]
__license__ = "MIT"
__date__ = "23/11/2016"

from setuptools import setup

setup(name='dynXRD',
      version='0.0.1',
      description='Python dynamic crystal diffraction calcution',
      author='Federica Maria Surace, Carsten Richter, Dmitri Novikov, Manuel Sanchez del Rio',
      author_email='federicamsurace@gmail.com',
      url='https://github.com/suracefm/dynXRD/',
      packages=['dynXRD',
                'dynXRD.tests',],
      install_requires=[
                        'numpy',
                        'sympy',
                        #'xraylib',
                       ],
      # test_suite='tests'
     )

