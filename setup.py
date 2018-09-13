#!/usr/bin/env python
from setuptools import setup

version = '0.1'
name = 'neoantigens'

setup(
    name=name,
    version=version,
    author='Vlad Saveliev',
    description='UMCCR neoantigens calling pipeline',
    keywords='bioinformatics',
    license='GPLv3',
    packages=[
        name,
    ],
    include_package_data=True,
)
