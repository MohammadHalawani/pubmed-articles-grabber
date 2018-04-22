from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pubmed-articles-grabber',  # Required
    version='0.0.1',  # Required
    description='Grabs and downloads full-text versions of PubMed records in different formats',  # Required

    long_description=long_description,  # Optional

    long_description_content_type='text/markdown',  # Optional (see note above)

    url='https://github.com/ezzohala/pubmed-articles-grabber',  # Optional

    author='Mohammad K. Halawani',  # Optional

    author_email='mkhalawani@uqu.edu.sa',  # Optional

    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
    ],

    keywords='pubmed grabber full-text article-extractor pubmed-central doi crossref-api',  # Optional

    packages=find_packages(exclude=['contrib', 'docs', 'tests']),  # Required

    install_requires=['unicodecsv', 'pandas', 'requests', 'habanero', 'biopython'],  # Optional


    # If there are data files included in your packages that need to be
    # installed, specify them here.
    #
    # If using Python 2.6 or earlier, then these have to be included in
    # MANIFEST.in as well.
    package_data={  # Optional
        'clickThroughToken': ['clickThroughToken.txt'],
        'wantedListOfPMIDs': ['wanted.csv'],
    },

    project_urls={  # Optional
        'Source': 'https://github.com/ezzohala/pubmed-articles-grabber',
    },
)

