#!/usr/bin/env python

from argparse import ArgumentParser
import cProfile
import logging
import numpy as np
import os
import sys


def script(data_path, out_path, num_diseases, **kwargs):
    cwd = os.getcwd()
    try:
        os.makedirs(out_path)
    except OSError:
        if not os.path.isdir(out_path):
            raise

    mean_reciprocal_rank_dict = {}
    total_count = 0

    # Assume the results data path contains folders for different experimental
    # conditions
    for experimental_condition in os.listdir(data_path):
        rank_counts = np.zeros(num_diseases)
        for patient_result_file in os.listdir(os.path.join(data_path,
                                                           experimental_condition)):
            name, extension = os.path.splitext(patient_result_file)
            if extension == '.rank':
                with open(os.path.join(cwd, data_path, experimental_condition,
                                       patient_result_file)) as f:
                    rank = int(list(f)[0].split('\t')[0])
                    try:
                        rank_counts[rank] += 1
                        total_count += 0
                    except IndexError:
                        print 'Ignored a rank count because it exceeded the'
                        'specified rank consideration range:'
                        '\t' + experimental_condition + ':', rank

        # mean reciprocal rank = harmonic mean of ranks
        mean_reciprocal_rank = np.sum(np.reciprocal([r for r in rank_counts
                                                        if r > 0]))
        mean_reciprocal_rank = np.divide(mean_reciprocal_rank, total_count)
        mean_reciprocal_rank_dict[experimental_condition] = str(mean_reciprocal_rank)

        # ROC and precision / recall stuff
        # 'threshold rates' are the ranks
        # so the (FP, TP) = (0, 0) classifier is the rank 0 threshold classifier
        # and the (FP, TP) = (1, 1) classifier is the rank infty threshold classifier

        # create the lines ahead of time so we can get rid of duplicates
        ROC_to_write = []
        PR_to_write = []

        true_pos = 0
        false_neg = np.sum(rank_counts)
        false_pos = 0
        true_neg = np.sum([j * rank_counts[j] for j in range(len(rank_counts))])

        for i in range(len(rank_counts)):

            # brute force computation
            #true_pos = np.sum(rank_counts[:i])
            #false_neg = np.sum(rank_counts[i:])
            #false_pos = np.sum([(j-1)*rank_counts[j] for j in
            #                    range(len(rank_counts[:i]))])
            #true_neg = np.sum([(j)*rank_counts[j] for j in
            #                   range(i, len(rank_counts))])

            true_pos += rank_counts[i-1]
            false_neg -= rank_counts[i-1]
            false_pos += (i-2) * rank_counts[i-1]
            true_neg -= (i-2) * rank_counts[i-1]

            recall = true_pos / (true_pos + false_neg) # true positive rate
            precision = true_pos / (true_pos + false_pos) # positive pred value
            specificity = true_neg / (false_pos + true_neg) # true negative rate
            fall_out = false_pos / (false_pos + true_neg) # false positive rate

            ROC_to_write.append((str(fall_out) + ',' +  str(recall) + '\n'))
            PR_to_write.append((str(recall) + ',' +  str(precision) + '\n'))

        ROC_to_write = list(set(ROC_to_write))
        ROC_to_write.sort()
        PR_to_write = list(set(PR_to_write))
        PR_to_write.sort()

        name = experimental_condition + '_ROC.analysis'
        with open(os.path.join(cwd, out_path, name), 'w') as f:
            f.write('fall-out,recall\n')
            f.writelines(ROC_to_write)

        name = experimental_condition + '_PR.analysis'
        with open(os.path.join(cwd, out_path, name), 'w') as f:
            f.write('recall,precision\n')
            f.writelines(PR_to_write)

    with open(os.path.join(cwd, out_path, 'mean_reciprocal_rank.analysis'), 'w') as f:
        f.write('condition,mean reciprocal rank\n')
        for experimental_condition in os.listdir(data_path):
            f.write(','.join([experimental_condition,
                              mean_reciprocal_rank_dict[experimental_condition] + '\n']))

def parse_args(args):
    parser = ArgumentParser()

    parser.add_argument('--data_path', '-D', metavar='DATA', required=True,
            help='Directory from which to grab results data files')
    parser.add_argument('--out_path', '-O', metavar='OUT', required=True,
            help='Directory in which to write results')
    parser.add_argument('--num_diseases', '-N', metavar='NUM', type=int,
            default=500, help='Number of diseases used to generate results')
    parser.add_argument('--logging', default='WARNING',
            choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
            help='logging level')
    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.logging)
    script(**vars(args))


if __name__ == '__main__':
    cProfile.run('sys.exit(main())')
