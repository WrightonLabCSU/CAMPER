from setuptools import setup, find_packages
from camper_dramkit import __version__ as version
import os

__author__ = 'rmflynn'
__version__ = version


here = os.path.abspath(os.path.dirname(__file__))
data_folder = os.path.join(here, 'camper_dramkit', 'data')
if not os.path.exists(data_folder):
    os.mkdir(data_folder)
    os.rename(os.path.join(here, '..', 'CAMPER_blast.faa'),
              os.path.join(data_folder, 'CAMPER_blast.faa'))
    os.rename(os.path.join(here, '..', 'CAMPER_blast_scores.tsv'),
              os.path.join(data_folder, 'CAMPER_blast_scores.tsv'))
    os.rename(os.path.join(here, '..', 'CAMPER.hmm'),
              os.path.join(data_folder, 'CAMPER.hmm'))
    os.rename(os.path.join(here, '..', 'CAMPER_distillate.tsv'),
              os.path.join(data_folder, 'CAMPER_distillate.tsv'))
    os.rename(os.path.join(here, '..', 'CAMPER_hmm_scores.tsv'),
              os.path.join(data_folder, 'CAMPER_hmm_scores.tsv'))

with open(os.path.join(here, "..", 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="camper_dramkit",
    version=__version__,
    description="A tool to use the CAMPER dataset",
    long_description=long_description,
    long_description_content_type='text/markdown',  # Optional (see note above)
    packages=['camper_dramkit'],
    package_dir={'camper_dramkit': 'camper_dramkit'},
    package_data={'camper_dramkit': ['data/*']},
    python_requires='>=3.8',
    install_requires=['pandas', 'altair', 'sqlalchemy', 'networkx', 'openpyxl', 'numpy', 'click'],
    entry_points={
        'console_scripts': [
            'camper_distill = camper_dramkit.camper_distill:summarize_genomes_cmd',
            'camper_annotate = camper_dramkit.camper_annotate:annotate_genes_cmd',
            'combine_annotations_lowmem = camper_dramkit.combine_annotations_lowmem:append_annotations_lowmem'
        ],
    },
    author="Rory Flynn",
    author_email='rory.flynn@colostate.edu',
    url="https://github.com/WrightonLabCSU/CAMPER",
    download_url="https://github.com/WrightonLabCSU/CAMPER/tarball/%s" % __version__
)
