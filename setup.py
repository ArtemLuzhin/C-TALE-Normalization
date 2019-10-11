#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

import re
VERSIONFILE="ctale_normalize/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(
      name='CTALE_normalize',
      version=verstr,
      packages=['ctale_normalize'],
      entry_points={
          'console_scripts': ['ctale_normalize = ctale_normalize.__main__:main']},
      install_requires=['Cython', 'numpy', 'cooler', 'pandas', 'natsort',
                        'scipy', 'cooltools'],
      description='Normalize your C-TALE data!',
      long_description=long_description,
      long_description_content_type='text/markdown',
      project_urls={'Source':'https://github.com/ArtemLuzhin/C-TALE-Normalization',
                    'Issues':'https://github.com/ArtemLuzhin/C-TALE-Normalization/issues'},
      author='Artem Luzhin',
      author_email='artyom.luzhin@gmail.com',
      classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)
