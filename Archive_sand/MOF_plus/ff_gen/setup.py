#!/usr/bin/env python
# -*- coding: utf-8 -*-

from glob import glob
import os, sys
from distutils.core import setup
from distutils.command.install_data import install_data

setup(
    name='ff_gen',
    version='2.0.1',
    description='Python library to derive MOF-FF force fields',
    author='Johannes P. Dürholt',
    author_email='johannes.duerholt@rub.de',
    package_dir = {'ff_gen': 'ff_gen'},
    packages=['ff_gen'],
    scripts=['scripts/atomtyper',
            'scripts/check_closefit',
            'scripts/collect_data.py',
            'scripts/create_key',
            'scripts/create_ref',
            'scripts/diffkey',
            'scripts/initial',
            'scripts/txyz2xyz',
            'scripts/xyz2txyz'
            ]
)
