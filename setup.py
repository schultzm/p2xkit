from setuptools import setup, find_packages
import p2xkit
import os

def read(fname):
    '''
    Read the README
    '''
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'p2xkit',
    version = p2xkit.__version__,
    description = p2xkit.__description__,
    long_description=read('README.md'),
    classifiers = ['Development Status :: 3 - Alpha',
                   'License :: OSI Approved :: GNU Affero General ' +
                   'Public License v3 or later (AGPLv3+)',
                   'Programming Language :: Python :: 3.5',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Medical Science Apps.',
                   'Intended Audience :: Science/Research',
                   "Programming Language :: Python :: 3.7"],
    keywords = ["primer", "oligo"],
    download_url = p2xkit.__download_url__,
    author = p2xkit.__author__,
    author_email = p2xkit.__author_email__,
    license = p2xkit.__license__,
    packages=find_packages(exclude=["contrib", "docs", "tests"]),  # Required
    python_requires=">=3.7",
    include_package_data = True,
    install_requires = ["biopython>=1.76",
                         "pysam>=0.15.4"],
    extras_require={  # Optional
        "dev": ["pre-commit", "pipenv"],
        "test": ["pytest", "pytest-cov"],
    },
    package_data={"": [""]},  # Optional
    entry_points={"console_scripts": ["p2xkit = p2xkit.__main__:main"]}  # Optional

    )
