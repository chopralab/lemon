from skbuild import setup

LONG_DESCRIPTION = """Lemon is a framework and API for mining information from
the PDB and other sources. This module provides the Python bindings to the Lemon
framework to allow users to create Lemon workflows completly in Python without
building any C++ executable. See the documentation at
http://chopralab.github.io/lemon"""

setup(
    name= "lemon_proteins",
    version="0.0.1",
    long_description=LONG_DESCRIPTION,
    description="Mine data from the PDB in minutes",
    keywords="chemistry computational cheminformatics proteins structural biology",
    author="Jonathan Fine",
    author_email="choprait@purdue.edu",
    license="BSD",
    url="http://github.com/chopralab/lemon",
    setup_requires=["scikit-build"],
    cmake_args=[
        '-DLEMON_BUILD_PYTHON:BOOL=OFF',
        '-DLEMON_LINK_SHARED:BOOL=OFF',
        '-DLEMON_BUILD_PROGS:BOOL=OFF',
    ],
    packages=['lemon_proteins'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: Unix",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
