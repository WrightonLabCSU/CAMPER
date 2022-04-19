from setuptools import setup, find_packages
from mag_annotator import __version__ as version
from os import path

__author__ = 'shafferm'
__version__ = version

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="DRAM-bio",
    version=__version__,
    # scripts=['scripts/DRAM.py', 'scripts/DRAM-v.py', 'scripts/DRAM-setup.py'],
    packages=find_packages(),
    description="A tool to use the CAMPER dataset",
    long_description=long_description,
    long_description_content_type='text/markdown',  # Optional (see note above)
    package_data={'mag_annotator': ['CONFIG']},
    python_requires='>=3.8',
    install_requires=['scikit-bio', 'pandas', 'altair', 'sqlalchemy', 'networkx', 'openpyxl', 'numpy', 'click'],
    entry_points={
        'console_scripts': [
            'camper_distill = camper_distill:summarize_genomes',
            'camper_annotate = camper_annotate:annotate_genes',
            'combine_annotations_lowmem = combine_annotations_lowmem:append_annotations_lowmem',
        ],
    },
    author="Michael Shaffer",
    author_email='michael.t.shaffer@colostate.edu',
    url="https://github.com/shafferm/DRAM/",
    download_url="https://github.com/shafferm/DRAM/tarball/%s" % __version__
)
