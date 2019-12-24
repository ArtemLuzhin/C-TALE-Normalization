#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ctale_normalize import CTALE_norm_iterative, save_weight
import argparse
import cooler
import logging
from natsort import natsorted

def main():
    parser = argparse.ArgumentParser(
                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("cooler", type=str,
                        help="Cooler file with your C-TALE data")
    parser.add_argument("coordinates", type=str,
                        help="""Coordinates of the captured region in UCSC
                        format. If multiple regions (one per chromosome!) were
                        captured, semicolon-separated coordinates for each of
                        them can be used. For example,
                        chr1:150,000-151,000;chr2:320000-420000""")
    parser.add_argument("--mult_factor", type=float, default=1.54,
                        required=False,
                        help="Factor for correction of zone 3")
    parser.add_argument("--IC_steps", type=int, default=50,
                        required=False,
                        help="""Number of steps for iterative correction in
                                zone 2""")
    parser.add_argument("--tolerance", type=float, default=10**-5,
                        required=False,
                        help="""Target variance for iterative correction""")
    parser.add_argument("--MAD_max", type=float, default=4,
                    required=False,
                    help="""Median absolute deviation filter value to remove
                            very low/high covered bins""")
#    parser.add_argument("output", type=str,
#                        help="Where to save the output")

    args = parser.parse_args()

    logging.basicConfig(format='%(message)s',
                        level=getattr(logging, 'INFO'))

    C=cooler.Cooler(args.cooler) #load coolfile
    logging.info('Loaded cool')
    regions = args.coordinates.split(';')
    chroms = []
    weights = []
    bal_info = []
    for region in natsorted(regions):
        chrom, start, end = cooler.util.parse_region(region)
        if chrom in chroms:
            raise ValueError("""Two regions are on the same chromosome, this is
                                not yet supported""")
        #load raw matrix
        mtx = C.matrix(balance=False, sparse=True).fetch(chrom).tocsr()
        logging.info('Loaded matrix for %s' % chrom)
        #Perform normalization
        logging.info('Normalization...')
        weight, converged = CTALE_norm_iterative(mtx, start, end,
                                                        C.binsize,
                                                        steps=args.IC_steps,
                                                        mult=args.mult_factor,
                                                        tolerance=args.tolerance,
                                                        mad_cutoff=args.MAD_max)
        info = {'tol':args.tolerance,
                'mad_max':args.MAD_max,
                'converged':converged}
        chroms.append(chrom)
        weights.append(weight)
        bal_info.append(info)
    save_weight(args.cooler, chroms, weights, bal_info)
#    logging.info('Saving as '+ args.output)
    #Save_coolfile
#    logging.info('Finished!')