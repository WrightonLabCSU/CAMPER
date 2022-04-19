from setuptools import setup, find_packages
from camperkit import __version__ as version
from os import path

__author__ = 'rmflynn'
__version__ = version

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="CAMPER_DRAMKit",
    version=__version__,
    description="A tool to use the CAMPER dataset",
    long_description=long_description,
    long_description_content_type='text/markdown',  # Optional (see note above)
    packages=['camperkit'],
    package_data={'camperdb': ['../CAMPERdb/*']},
    python_requires='>=3.8',
    install_requires=['scikit-bio', 'pandas', 'altair', 'sqlalchemy', 'networkx', 'openpyxl', 'numpy', 'click'],
    entry_points={
        'console_scripts': [
            'camper_distill = camperkit.camper_distill:summarize_genomes',
            'camper_annotate = camperkit.camper_annotate:annotate_genes',
            'combine_annotations_lowmem = camperkit.combine_annotations_lowmem:append_annotations_lowmem'
        ],
    },
    author="Rory Flynn",
    author_email='rory.flynn@colostate.edu',
    url="https://github.com/WrightonLabCSU/CAMPER",
    download_url="https://github.com/WrightonLabCSU/CAMPER/tarball/%s" % __version__
)
