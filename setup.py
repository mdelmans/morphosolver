from setuptools import setup, find_packages
from codecs import open
from os import path

setup(
	name = 'morphosolver',
	version = '1.0',
	description='morphosolver',
	author = 'Mihails Delmans',
	author_email='md656@cam.ac.uk',
	license = 'MIT',
	classifiers=[
		'Development Status :: 4 - Beta',
		'Intended Audience :: Developers',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2.7',
	],
	packages = find_packages()
)

