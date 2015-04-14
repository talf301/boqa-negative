#!/usr/bin/env python

import logging
import sys
import os
from argparse import ArgumentParser

from net import Net

__author__ = 'Tal Friedman (talf301@gmail.com)'

def script(data_path, patient_path, out_path, p_sampling, sampling, k_freqs, **kwargs):
    try:
        net = Net(os.path.join(data_path, '/hp.obo'), os.path.join(data_path, '/phenotype_annotations.tab'), os.path.join(data_path, '/negative_phenotype_annotation.tab'))
    except IOError, e:
        logging.error(e)
        sys.exit(1)

    # Figure out if we have a single patient or a directory of them
    hpo_files = []
    try:
        files = os.listdir(patient_path)
        hpo_files = [f for f in files if f.endswith('_hpo.txt')]
    except OSError:
        hpo_files.append(patient_path)

    #for hpo_file in hpo_files:


def parse_args(args):
    parser = ArgumentParser()

    parser.add_argument('--data_path', '-D', metavar='DATA', default='./',
            help='Directory from which to grab data files')
    parser.add_argument('--patient_path', '-P', metavar='PATIENTS', required=True,
            help='Directory from which to grab patients in _hpo.txt format, or a single'
            "_hpo.txt file")
    parser.add_argument('--out_path', '-O', metavar='OUT', required=True,
            help='Directory in which to dump results files from diagnosis')
    parser.add_argument('--logging', default='WARNING',
            choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
            help='logging level')
    parser.add_argument('--p_sampling', '-p', metavar='P_SAMPLING', type=float, default=0,
            help='Do inference via p sampling, giving the value of p. This will make'
                 'use of negative phenotype annotations')
    parser.add_argument('--sampling', '-s', metavar='SAMPLING', type=int, default=0,
            help='Do inference via sampling, indicate the number of samples as an integer'
            'greater than 0.')
    parser.add_argument('--k_freqs', '-k', metavar='K_FREQ', type=int, default=0,
            help='Do inference by looking at the K lowest frequency annotations and'
            'expanding by expectation')
    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.logging)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
