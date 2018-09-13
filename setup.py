#!/usr/bin/env python
from setuptools import setup
from os.path import join

version = '0.2'
name = 'nag'

setup(
    name=name,
    version=version,
    author='Vlad Saveliev',
    description='UMCCR neoantigens calling pipeline',
    keywords='bioinformatics',
    license='GPLv3',
    packages=[name],
    scripts=[
        join(name, 'pizzly_to_bedpe.py'),
        join(name, name),
    ],
    include_package_data=True,
)
