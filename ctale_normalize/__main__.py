#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ctale_normalize import CTALE_norm_iterative, Save_coolfile
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
                        captured, write comma-separated coordinates for each
                        of them""")
    parser.add_argument("--mult_factor", type=float, default=1.54,
                        required=False,
                        help="Factor for correction of zone 3")
    parser.add_argument("--IC_steps", type=int, default=20,
                        required=False,
                        help="""Number of steps for iterative correction in
                                zone 2""")
    parser.add_argument("output", type=str,
                        help="Where to save the output")

    args = parser.parse_args()

#    logging.basicConfig(format='%(message)s',
#                        level=getattr(logging, args.logLevel))

    C=cooler.Cooler(args.cooler) #load coolfile
    logging.info('Loaded cool')
    regions = args.coordinates .split(',')
    chroms = []
    mtxs = []
    for region in natsorted(regions):
        chrom, start, end = cooler.util.parse_region(region)
        if chrom in chroms:
            raise ValueError("""Two regions are on the same chromosome, this is
                                not yet supported""")
        #load raw matrix
        mtx=C.matrix(balance=False, sparse=True).fetch(chrom)
        logging.info('Loaded matrix for %s' % chrom)
        #Perform normalization
        logging.info('Normalization...')
        mtxs.append(CTALE_norm_iterative(mtx, start, end, C.binsize, steps=args.IC_steps,
                               mult=args.mult_factor))
        chroms.append(chrom)
    logging.info('Save as '+ args.output)
    #Save_coolfile
    Save_coolfile(C, chroms, mtxs, args.output)
    logging.info('Finished!')