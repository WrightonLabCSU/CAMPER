from setuptools import setup, find_packages
from camper_dramkit import __version__ as version
from os import path

__author__ = 'rmflynn'
__version__ = version

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="camper_dramkit",
    version=__version__,
    description="A tool to use the CAMPER dataset",
    long_description=long_description,
    long_description_content_type='text/markdown',  # Optional (see note above)
    packages=['camper_dramkit'],
    include_package_data=True,
    # data_files=[
    #     ('camper_dramkit/CAMPERdb', ["../CAMPER.hmm", "../CAMPER_blast_scores.tsv", 
    #                   "../CAMPER_distillate.tsv", "../CAMPER_hmm_scores.tsv"])
    # ],
    package_data={'camper_dramkit': ["../CAMPER.hmm", "../CAMPER_blast_scores.tsv", 
                      "../CAMPER_distillate.tsv", "../CAMPER_hmm_scores.tsv"]},
    python_requires='>=3.8',
    # install_requires=['pandas', 'altair', 'sqlalchemy', 'networkx', 'openpyxl', 'numpy', 'click', 'scikit-bio'],
    install_requires=['pandas', 'altair', 'sqlalchemy', 'networkx', 'openpyxl', 'numpy', 'click', 'scikit-bio'],
    entry_points={
        'console_scripts': [
            'camper_distill = camper_dramkit.camper_distill:summarize_genomes',
            'camper_annotate = camper_dramkit.camper_annotate:annotate_genes',
            'combine_annotations_lowmem = camper_dramkit.combine_annotations_lowmem:append_annotations_lowmem'
        ],
    },
    author="Rory Flynn",
    author_email='rory.flynn@colostate.edu',
    url="https://github.com/WrightonLabCSU/CAMPER",
    download_url="https://github.com/WrightonLabCSU/CAMPER/tarball/%s" % __version__
)
