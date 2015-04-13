#!/usr/bin/env python


import logging
import os
import random
import sys

import hpo

from argparse import ArgumentParser
from omim import MIM


__author__ = 'Tal Friedman (talf301@gmail.com)'


def add_noise(num, phenotypes, omim_dict):
    """
    Randomly add terms to the given list of phenotypes.

    Args:
        num: number of random phenotypes to add
        phenotypes: list of phenotypes to add noise to
        omim_dct: a dict of OMIM number -> omim.Disease

    Returns:
        A noisy list of phenotypes.

    """
    # List of phenotypes to sample from
    sample_list = [x for y in omim_dict.itervalues() for x in y.phenotype_freqs]

    # Do sampling
    sampled_pheno = random.sample(sample_list, num)
    # Finally combine them and remove repetition
    phenotypes.extend(sampled_pheno)
    phenotypes = list(set(phenotypes))

    return phenotypes

def add_imprecision(hp, phenotypes):
    """
    Randomly add some imprecision by moving terms up the ontology tree at random
    to simulate different language usage.

    Args:
        hp: an hpo.HPO instance
        phenotypes: list of phenotypes to add imprecision to

    Returns:
        A list of the phenotypes with imprecision added.
        Note that the list may be shorter if intersection occurs.

    """
    # We're going to use a set, since scaling phenotypes up could
    # potentially make them intersect
    new_pheno = set()

    for pheno in phenotypes:
        try:
            # Key may not be found since there are things like
            # inheritance patterns in phenotypic annotations
            ancestors = list(hpo.get_ancestors(hp[pheno]))
        except KeyError:
            continue

        # Get as string ids and remove root
        ancestors = map(lambda x: x.id, ancestors)
        ancestors.remove('HP:0000118')

        # Now randomly add an ancestor
        new_pheno.add(random.choice(ancestors))

    # Return phenotypes as a list
    return list(new_pheno)

def sample_phenotypes(pos_omim_dict, neg_omim_dict, disease, hp, imprecision,
                      noise, default_freq=1.0):
    """
    Sample phenotypes randomly from a disease.

    Args:
        pos_omim_dict: dict of OMIM number -> omim.Disease containing positive
        annotations
        neg_omim_dict: dict of OMIM number -> omim.Disease containing negative
        annotations
        disease: an OMIM number representing the disease
        hp: an hpo.HPO instance
        imprecision: whether or not to add imprecision
        default_freq: default frequency to use if not specified

    Each disease should have at least one phenotype entry.

    Returns:
        A list of sampled phenotypes.

    """
    phenotypes = []

    # Grab the positive omim.Disease for this OMIM ID number
    try:
        pos_omim_dis = pos_omim_dict[disease]
    except KeyError:
        logging.warning('Could not find OMIM entry for %s' % disease)

    # Grab the negative omim.Disease for this OMIM ID number
    try:
        neg_omim_dis = neg_omim_dict[disease]
    except KeyError:
        logging.warning('Could not find negative annotations for ID %s' % disease)

    # If frequency available, we will sample, otherwise we use default freq
    phenotype_freqs = pos_omim_dis.phenotype_freqs
    assert phenotype_freqs, "Missing phenotypes for: %s" % disease
    for pheno, freq in phenotype_freqs.iteritems():
        if not freq and random.random() < default_freq:
            phenotypes.append(pheno)
        else:
            if random.random() < freq:
                phenotypes.append(pheno)
    if phenotypes:
        # Log the original number of phenotypes
        orig_len = len(phenotypes)
        # Add imprecision if necessary
        if imprecision:
            phenotypes = add_imprecision(hp, phenotypes)
        # Add noise if necessary
        if noise:
            phenotypes = add_noise(int(orig_len * noise), phenotypes, pos_omim_dict)
        return phenotypes
    else:
        logging.warning("Random phenotype sampling for %s resulted in"
                " empty set" % disease)
        return sample_phenotypes(pos_omim_dict, neg_omim_dict, disease, hp, imprecision, noise,
                                 default_freq)

def infect_pheno(patient, disease, pos_omim_dict, neg_omim_dict, hp, imprecision, noise, default_freq):
    """
    Do phenotypic infection by producing a file with a list of phenotypes
    associated with disease.

    Args:
        patient: a string path to the patient vcf
        disease: an OMIM disease ID to infect patient
        pos_omim_dict: a dict of OMIM number -> omim.Disease
        pos_omim_dict: a dict of OMIM number -> omim.Disease
        hp: an hpo.HPO instance
        imprecision: whether or not to add imprecision to pheno sampling
        default_freq: default frequency to use if info not found

    """
    # Sample phenotypes
    phenotypes = sample_phenotypes(pos_omim_dict, neg_omim_dict, disease, hp, imprecision, noise,
                                   default_freq)

    assert patient.endswith('.vcf')

    # If patient is HG01.vcf, phenotypes in HG01_hpo.txt
    pheno_file = patient[:-4] + '_hpo.txt'
    with open(pheno_file, 'w') as hpo:
        hpo.write(','.join(phenotypes))

def load_data(data_path):
    """
    Load all the required data files the program needs.

    Args:
        data_path: string file path to the directory files are in

    Returns
        (OMIM number -> omim.Disease, OMIM number -> omim.Disease, HPO,
        list of str)

    """
    # Load hpo
    hp = hpo.HPO(os.path.join(data_path, 'hp.obo'))
    # Filter to phenotypic abnormailities
    hp.filter_to_descendants('HP:0000118')

    # Load and filter OMIM into an OMIM number: Disease dict,
    # Drop everything that isn't a child of 118
    pos_mim = MIM(os.path.join(data_path, 'phenotype_annotation.tab'))
    pos_omim = filter(lambda d:d.db == 'OMIM', pos_mim.diseases)
    for o in pos_omim:
        o.phenotype_freqs = {pheno:freq for pheno,freq in o.phenotype_freqs.iteritems() if hp.hps.has_key(pheno)}
    pos_omim_dict = {dis.id:dis for dis in pos_omim}

    # Load the negative annotations
    neg_mim = MIM(os.path.join(data_path, 'negative_phenotype_annotation.tab'))
    neg_omim = filter(lambda d:d.db == 'OMIM', neg_mim.diseases)
    for o in neg_omim:
        o.phenotype_freqs = {pheno:freq for pheno,freq in o.phenotype_freqs.iteritems() if hp.hps.has_key(pheno)}
    neg_omim_dict = {dis.id:dis for dis in neg_omim}

    # Get the disease possibilities
    diseases = [d.id for d in pos_mim.diseases if d.db == 'OMIM']

    # There is one disease with only negative annotations
    # Currently ignored
    #assert set([d.id for d in neg_mim.diseases if d.db == 'OMIM']).issubset(diseases)

    return pos_omim_dict, neg_omim_dict, hp, diseases

def script(data_path, out_path, generate, num_samples, default_freq,
        imprecision, noise, negative_phenotypes, **kwargs):
    try:
        pos_omim_dict, neg_omim_dict, hp, diseases = load_data(data_path)
    except IOError, e:
        logging.error(e)
        sys.exit(1)

    # Dealing with pairs
    if generate == 'PAIRS':
        for i in range(num_samples):
            # First, get a disease
            # A uniform sample over all diseases
            if negative_phenotypes:
                found = False
                while not found:
                    disease = random.choice(diseases)
                    try:
                        neg_omim_dict[disease]
                        found = True
                    except KeyError:
                        continue
            else:
                disease = random.choice(diseases)

            # Name patients based on disease and iteration
            new_pair = ['First_' + disease + '_' + str(i) + '.vcf',
                    'Second_' + disease + '_' + str(i) + '.vcf']

            # Finally, infect both patients with disease
            for patient in new_pair:
                infect_pheno(os.path.join(out_path, patient), disease, pos_omim_dict,
                    neg_omim_dict, hp, imprecision, noise, default_freq)

    # Dealing with individual patients
    if generate == 'PATIENTS':
        for i in range(num_samples):
            # First, get a disease
            # A uniform sample over all diseases
            disease = random.choice(diseases)

            # Name patient based on disease and iteration
            new_patient = disease + '_' + str(i) + '.vcf'

            # Finally, infect patient with disease
            infect_pheno(os.path.join(out_path, new_patient), disease, pos_omim_dict,
                neg_omim_dict, hp, imprecision, noise, default_freq,
                negative_phenotypes)

def parse_args(args):
    parser = ArgumentParser()

    parser.add_argument('--data_path', '-d', metavar='DATA', default='./',
            help='Directory from which to grab data files')
    parser.add_argument('--out_path', '-o', metavar='OUT', required=True,
            help='Directory in which to put the generated patient files')
    parser.add_argument('--generate', dest='generate',
            choices=['PATIENTS', 'PAIRS'], default='PAIRS',
            help='Generate pairs or individuals (default is pairs)')
    parser.add_argument('-N', type=int, dest='num_samples', required=True,
            help='Number of samples to generate (patients or pairs)')
    parser.add_argument('-D', '--default_freq', type=float,
            help='Default frequency for phenotype if info not found (default'
            'is 1.0)', default=1.0)
    parser.add_argument('--imprecision', action='store_true',
            help='Add imprecision of hpo term selection (randomly send'
            'terms to ancestors')
    parser.add_argument('--noise', type=float, default=0.0,
            help='Add phenotypic noise (random phenotypes). Recommended'
            'amount is 0.5 (i.e., half of the real amount will be added'
            'as noise).')
    parser.add_argument('--negative_phenotypes', action='store_true',
            help='Restrict sampled diseases to those that have negative'
            'annotations.')
    parser.add_argument('--logging', default='WARNING',
            choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
            help='Logging level')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.logging)
    script(**vars(args))


if __name__ == '__main__':
    sys.exit(main())
