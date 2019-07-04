#!/usr/bin/env python
# coding: utf-8
import os
import papermill as pm
import sys
import subprocess
import argparse


def main():
    # parse_args
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample_list', nargs='+',
                        help='<Required> sample list', required=True)
    parser.add_argument('-n', '--sample_name_list', nargs='+',
                        help='<Required> sample name list', required=True)
    parser.add_argument('-o', '--output_dir',
                        help='<Required> output dir', required=True)
    all_args = parser.parse_args()

    sample_ids = all_args.sample_list
    sample_names = all_args.sample_name_list
    output_dir = all_args.output_dir
    output_file = os.path.join(output_dir, 'clustering_merge.ipynb')

    # excute  jupyter notebook
    print("Running jupyter:", sample_names)
    pm.execute_notebook(
        'clustering_merge.ipynb',
        output_file,
        parameters={"samples": sample_ids,
                    "sample_names": sample_names,
                    }
    )

    subprocess.call(["jupyter-nbconvert", "--to", "html",
                     output_file])


if __name__ == '__main__':
    main()
