from __future__ import division
import pandas as pd
import numpy as np
import cooler
import scipy.stats
from scipy import sparse, stats
import sys
from natsort import natsorted
import logging

#functions

def CTALE_norm_multiplicate(mtx, ROI_start, ROI_end, resolution,
                            func=scipy.stats.gmean, mult=1.54):
    """mtx- matrix of individual chromosome/region +/- distance
    ROI_start - first coordinate of C-TALE region(bp)
    ROI_end - last coordinate of C-TALE region(bp)
    resolution - C-TALE map resolution(bp)
    func - function for calculating mean, can be numpy.mean,numpy.median and etc. Default: scipy.stats.gmean
    mult-coefficient of multiplication around ROI, default=1.54
    returns normalized matrix"""
    normalized=np.zeros(shape=mtx.shape)
    start_bin=int(ROI_start/resolution)
    end_bin=int(ROI_end/resolution)
    #fill main diagonal with zeros
    np.fill_diagonal(mtx,0)
    # first we multiply around region
    new_mtx=mtx*mult
    new_mtx[start_bin:end_bin,start_bin:end_bin]=mtx[start_bin:end_bin,start_bin:end_bin]
    for i in range(start_bin,end_bin):
        #left=new_mtx[i,i:]
        #up=new_mtx[:i,i]
        normalized[i,i:]=new_mtx[i,i:]/(np.sum(new_mtx[i,i:])+np.sum(new_mtx[:i,i])) #left
        normalized[:i,i]=new_mtx[:i,i]/(np.sum(new_mtx[:i,i])+np.sum(new_mtx[i,i:])) #up
    for i in range(start_bin,end_bin):
        for j in range(i+1,end_bin):
            normalized[i,j]=new_mtx[i,j]/func([np.sum(new_mtx[i,:]),np.sum(new_mtx[:,j])])

    i_lower = np.tril_indices(normalized.shape[0], -1) #creates symmetric matrix
    normalized[i_lower] = normalized.T[i_lower]
    return(normalized)

def multiplicate(mtx, ROI_start, ROI_end, resolution, mult=1.54):
    """mtx- matrix of individual chromosome/region +/- distance
    ROI_start - first coordinate of C-TALE region(bp)
    ROI_end - last coordinate of C-TALE region(bp)
    resolution - C-TALE map resolution(bp)
    mult-coefficient of multiplication around ROI, default=1.54
    Function perform multiplicating of matrix around ROI for selected coefficient"""
    start_bin = ROI_start//resolution
    end_bin = ROI_end//resolution
    #fill main diagonal with zeros
    mtx.setdiag(0)
    # first we multiply around region
    new_mtx = mtx*mult
    new_mtx[start_bin:end_bin+1, start_bin:end_bin+1] = mtx[start_bin:end_bin+1,
                                                            start_bin:end_bin+1]
    return new_mtx

def get_cov_var(mtx, start_bin, end_bin):
    cov1 = np.asarray(mtx.sum(axis=0)).ravel()
    cov2 = np.asarray(mtx.sum(axis=1)).ravel()

    cov1[start_bin:end_bin+1] /= cov1[start_bin:end_bin+1].mean()
    cov2[start_bin:end_bin+1] /= cov2[start_bin:end_bin+1].mean()

#    cov1mn = cov1[start_bin:end_bin+1].mean()
#    cov2mn = cov2[start_bin:end_bin+1].mean()

    cov1[:start_bin] = 1 #cov1mn
    cov1[end_bin+1:] = 1 #cov1mn

    cov2[:start_bin] = 1 #cov2mn
    cov2[end_bin+1:] = 1 #cov2mn

    cov = cov1+cov2#stats.gmean([cov1, cov2], axis=0)
    cov /= cov.mean()
    var = np.var(cov[start_bin:end_bin+1])

    return cov, var

def CTALE_norm(mtx, ROI_start, ROI_end, resolution):
    """Single iteration of CTALE balancing
    mtx - matrix of individual chromosome/region +/- distance
    ROI_start - first coordinate of C-TALE region(bp)
    ROI_end - last coordinate of C-TALE region(bp)
    resolution - C-TALE map resolution(bp)
    returns normalized matrix and variance at this iteration"""
    start_bin = ROI_start//resolution
    end_bin = ROI_end//resolution
    cov, var = get_cov_var(mtx, start_bin, end_bin)
    cov = (cov-1)*0.8 + 1
    mtx = mtx.multiply(1/cov).multiply(1/cov[np.newaxis].T).tocsr()
    cov, var = get_cov_var(mtx, start_bin, end_bin)
    return mtx, var

def CTALE_norm_iterative(mtx, ROI_start, ROI_end, resolution, mult=1.54,
                         steps=20, tolerance=10**-5):
    """Main function that perform normalization until variance>tolerance
    mtx- matrix of individual chromosome/region +/- distance
    ROI_start - first coordinate of C-TALE region(bp)
    ROI_end - last coordinate of C-TALE region(bp)
    resolution - C-TALE map resolution(bp)
    mult-coefficient of multiplication around ROI, default=1.54
    steps-number of iterations, by default=20
    tolerance-when variance<tolerance algorithm stops.
    returns normalized matrix"""
    start_bin = ROI_start//resolution
    end_bin = ROI_end//resolution
    out = multiplicate(mtx=mtx, ROI_start=ROI_start, ROI_end=ROI_end,
                   resolution=resolution, mult=mult)
    for i in range(steps):
        out, var = CTALE_norm(out, ROI_start, ROI_end, resolution)
        logging.info('Iteration %s: var: %s' % (i, var))
        if var < tolerance:
            logging.info('Variance below %s' % tolerance)
            out /= out.sum(axis=1)[start_bin:end_bin+1].mean()
            # Ensure zones 2 and 3 are scaled identically
#            newsum = out[start_bin:end_bin+1, start_bin:end_bin+1].sum()
#            oldsum = mtx[start_bin:end_bin+1, start_bin:end_bin+1].sum()
#            factor = newsum/oldsum
#            out *= factor
#            out[start_bin:end_bin+1, start_bin:end_bin+1] /= factor
            return out
    raise ValueError('Too many interation without convergence')


def get_pixels(mtx, zero_id=0):
    mtx_upper_diag = sparse.triu(mtx, k=0)
    sc = mtx_upper_diag.tocoo(copy=False)
    pixels = pd.DataFrame({'bin1_id': sc.row+zero_id,
                           'bin2_id': sc.col+zero_id,
                           'count': sc.data})
    return pixels

def Save_coolfile(coolfile, chroms, mtxs, output_coolfile):
    """Function change raw HiC matrix of cool file to user selected (normalized) and write it to new file.
    Because function rewrite data of original cool, later you should load it with balance=False flag.
    coolfile - original HiC file
    mtx - matrix to write
    output_coolfile - name of new cool file
    genome - genome assembly id"""
    #create bins
    bins = [coolfile.bins().fetch(chrom)[:] for chrom in chroms]
    bins = pd.concat(bins).reset_index(drop=True)
    zerobins = dict(bins.groupby('chrom').apply(lambda x: x.index.min()))
    #Create sparse matrix
    pixels = pd.concat([get_pixels(mtx, zerobins[chrom]) for mtx, chrom in zip(mtxs, chroms)])
    cooler.create_cooler(output_coolfile, bins, pixels,
                     assembly=coolfile.info[u'genome-assembly'],
                     dtypes={'count':float})
    return