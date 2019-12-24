from __future__ import division
import pandas as pd
import numpy as np
import cooler
from cooler import util
import scipy.stats
from scipy import sparse, stats
import sys
from natsort import natsorted
import logging
import h5py

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
        normalized[i,i:]=new_mtx[i,i:]/(np.nansum(new_mtx[i,i:])+np.nansum(new_mtx[:i,i])) #left
        normalized[:i,i]=new_mtx[:i,i]/(np.nansum(new_mtx[:i,i])+np.nansum(new_mtx[i,i:])) #up
    for i in range(start_bin,end_bin):
        for j in range(i+1,end_bin):
            normalized[i,j]=new_mtx[i,j]/func([np.nansum(new_mtx[i,:]),np.nansum(new_mtx[:,j])])

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

def get_cov(mtx, start_bin, end_bin, norm=True):
    mtx[mtx!=mtx] = 0
    cov1 = np.asarray(mtx.sum(axis=0)).ravel()
    cov2 = np.asarray(mtx.sum(axis=1)).ravel()

    if norm:
        cov1[start_bin:end_bin+1] /= np.nanmean(cov1[start_bin:end_bin+1])
        cov2[start_bin:end_bin+1] /= np.nanmean(cov2[start_bin:end_bin+1])

    cov1[:start_bin] = 1
    cov1[end_bin+1:] = 1

    cov2[:start_bin] = 1
    cov2[end_bin+1:] = 1

    cov = cov1+cov2
    if norm:
        cov /= np.nanmean(cov)

    return cov

#def CTALE_norm(mtx, ROI_start, ROI_end, resolution, bad_bins):
#    """Single iteration of CTALE balancing
#    mtx - matrix of individual chromosome/region +/- distance
#    ROI_start - first coordinate of C-TALE region(bp)
#    ROI_end - last coordinate of C-TALE region(bp)
#    resolution - C-TALE map resolution(bp)
#    returns normalized matrix and variance at this iteration"""
#    start_bin = ROI_start//resolution
#    end_bin = ROI_end//resolution
#    cov = get_cov(mtx, start_bin, end_bin, bad_bins)
#    cov = (cov-1)*0.8 + 1
#    return cov

def mad_max(coverage, cutoff=5):
    logNzMarg = np.log(coverage[coverage > 0])
    med_logNzMarg = np.median(logNzMarg)
    dev_logNzMarg = util.mad(logNzMarg)
    cutoff = np.exp(med_logNzMarg - cutoff * dev_logNzMarg)
    return np.where(coverage < cutoff)

def CTALE_norm_iterative(mtx, ROI_start, ROI_end, resolution, mult=1.54,
                         steps=20, tolerance=10**-5, mad_cutoff=5):
    """Main function that perform normalization until variance>tolerance
    mtx - matrix of individual chromosome/region +/- distance
    ROI_start - first coordinate of C-TALE region(bp)
    ROI_end - last coordinate of C-TALE region(bp)
    resolution - C-TALE map resolution(bp)
    mult - coefficient of multiplication around ROI, default=1.54
    steps - number of iterations, by default=20
    tolerance - when variance<tolerance algorithm stops.
    mad_cutoff - MAD filter value

    returns final weihgts, and whether the balancing converged"""
    start_bin = ROI_start//resolution
    end_bin = ROI_end//resolution
    out = multiplicate(mtx=mtx, ROI_start=ROI_start, ROI_end=ROI_end,
                   resolution=resolution, mult=mult)
    cov = get_cov(out, start_bin, end_bin)
    bad_bins = mad_max(cov[start_bin:end_bin+1], mad_cutoff)+start_bin
    weights = np.ones_like(cov)
    weights[bad_bins] = np.nan

    for i in range(steps):
        cov = get_cov(out.multiply(weights).multiply(weights[np.newaxis].T).tocsr(),
                               start_bin, end_bin)
        cov[bad_bins] = np.nan
        var = np.var(cov[start_bin:end_bin+1][~(bad_bins-start_bin)])
        logging.info('Iteration %s: var: %s' % (i, var))
        if var >= tolerance:
            cov = (cov-1)*0.8 + 1
            weights = weights/cov
            continue
        else:
            logging.info('Variance below %s' % tolerance)
            cov = get_cov(out, start_bin, end_bin, norm=False)[start_bin:end_bin+1]
            weights /= cov.mean()*2#**0.5
            weights = weights**0.5
#            out = out.multiply(weights).multiply(weights[np.newaxis].T).tocsr()
            return weights, True
    logging.warn('Too many iterations without convergence')
    return weights, False


#def get_pixels(mtx, zero_id=0):
#    mtx_upper_diag = sparse.triu(mtx, k=0)
#    sc = mtx_upper_diag.tocoo(copy=False)
#    pixels = pd.DataFrame({'bin1_id': sc.row+zero_id,
#                           'bin2_id': sc.col+zero_id,
#                           'count': sc.data})
#    return pixels

def save_weight(coolfile, chroms, weights, bal_info={}):
    """
    coolfile - File to add weights too
    chroms - List of chromosomes to modify
    weights - List of weight vectors to write to corresponding chromosomes
    bal_info - List of dicts to add as attributes of the weight column"""
    cool_path, group_path = util.parse_cooler_uri(coolfile)
    c = cooler.Cooler(coolfile)
    with h5py.File(cool_path, 'r+') as h5:
        grp = h5[group_path]
        # add the bias column to the file
        h5opts = dict(compression='gzip', compression_opts=6)
        allweights = np.asarray([np.nan]*grp['bins']['chrom'].shape[0])
        if 'weight' not in grp['bins'].keys():
            grp['bins'].create_dataset('weight', data=allweights, **h5opts)
        for chrom, weight, info in zip(chroms, weights, bal_info):
            s, e = c.extent(chrom)
            grp['bins']['weight'][s:e] = weight
            grp['bins']['weight'].attrs.update(info)

#        grp['bins']['weight'][grp['bins']['chrom']==chrom]=weights
    return