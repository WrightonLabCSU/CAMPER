import os
from glob import glob
import pandas as pd
import click

input = 'output/all_annotations/*/annotations.tsv'
output = 'results/all_annotations.tsv'
def append_annotations_lowmem(paths:list, output:str):
    paths = glob(input)
    all_headers = pd.concat([pd.read_csv(i, nrows=0, sep='\t', index_col=0) for i in paths])
    all_headers.to_csv(output, sep='\t')
    for i in paths:
        all_headers.append(
            pd.read_csv(i, sep='\t', index_col=0)
        ).to_csv(output, mode='a', header=False, sep='\t')


