#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ctale_normalize import CTALE_norm_iterative, Save_coolfile
import argparse
import cooler
import logging

def main():
    parser = argparse.ArgumentParser(
                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("cooler", type=str,
                        help="Cooler file with your C-TALE data")
    parser.add_argument("coordinates", type=str,
                        help="Coordinates of the captured region in ucsc format")
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
    chrom, start, end = cooler.util.parse_region(args.coordinates)
    #load raw matrix
    mtx=C.matrix(balance=False).fetch(chrom)
    logging.info('Loaded matrix...')
    #Perform normalization
    logging.info('Normalization...')
    N=CTALE_norm_iterative(mtx,start, end, C.binsize, steps=args.IC_steps,
                           mult=args.mult_factor)
    logging.info('Save as '+ args.output)
    #Save_coolfile
    Save_coolfile(C, N, chrom, args.output, C.info[u'genome-assembly'])
    logging.info('Finished!')