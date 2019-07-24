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
    parser.add_argument('-r', '--resultion', default=1,
                        help='<default=1> clustering resolution', required=True)
    parser.add_argument('--correct_batch',
                        dest='correct_batch', action='store_true', help='Correct batch (default)')
    parser.add_argument('--no_correct_batch',
                        dest='correct_batch', action='store_false', help='Not correct batch')
    parser.set_defaults(correct_batch=True)
    parser.add_argument('--correct_batch_by_rep',
                        dest='correct_batch_by_rep', action='store_true', help='Correct batch on rep')
    parser.set_defaults(correct_batch_by_rep=False)

    all_args = parser.parse_args()

    sample_ids = all_args.sample_list
    sample_names = all_args.sample_name_list
    output_dir = all_args.output_dir
    correct_batch = all_args.correct_batch
    correct_batch_by_rep = all_args.correct_batch_by_rep
    resultion = all_args.resultion

    nb_name = 'clustering_merge.ipynb' if correct_batch else 'clustering_merge_nobatch_correct.ipynb'
    if correct_batch_by_rep:
        nb_name = 'clustering_merge_correct_rep.ipynb'

    output_file = os.path.join(output_dir, nb_name)

    # excute  jupyter notebook
    print("Running", nb_name, 'for', sample_names)

    pm.execute_notebook(
        nb_name,
        output_file,
        parameters={"samples": sample_ids,
                    "sample_names": sample_names,
                    "output_dir": output_dir,
                    "res": str(resultion),
                    }
    )

    subprocess.call(["jupyter-nbconvert", "--to", "html",
                     output_file])


if __name__ == '__main__':
    main()
