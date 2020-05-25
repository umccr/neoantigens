#!/usr/bin/env python
from os.path import join
from setuptools import setup
import versionpy
package_name = 'nag'
version = versionpy.get_version(package_name)

setup(
    name=package_name,
    version=version,
    author='Vlad Saveliev',
    description='UMCCR neoantigens calling pipeline',
    keywords='bioinformatics',
    license='GPLv3',
    packages=[package_name],
    scripts=[
        join('scripts', 'pizzly_to_bedpe.py'),
        join('scripts', 'nag'),
    ],
    include_package_data=True,
)
