import numpy as np
import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
import seaborn as sns
from brainspace.plotting import plot_hemispheres
from brainspace.utils.parcellation import map_to_labels
import hcp_utils as hcp
from brainsmash.mapgen.stats import nonparp
from brainsmash.mapgen.stats import pearsonr
from brainsmash.mapgen.stats import spearmanr
from brainspace.null_models import SampledSurrogateMaps
from scipy.sparse.csgraph import dijkstra
from brainspace.datasets import load_conte69
from scipy.interpolate import make_interp_spline
from scipy.io import loadmat
from matplotlib.colors import ListedColormap
import math
lh, rh = load_conte69()

def spin_spearman(x,y):
  # x, y should be one array 
  n_surrogate_datasets = 1000

  # Note: number samples must be greater than number neighbors
  num_samples = 25
  num_neighbors = 30
  distance= dijkstra(np.array(pd.read_csv('../data/LeftParcelGeodesicDistmat.txt',
                                          header=None, delimiter=' ')), directed=False)
  distance_idx_sorted = np.argsort(distance, axis=1)
  ssm = SampledSurrogateMaps(ns=num_samples, knn=num_samples, random_state=0,resample=True)
  ssm.fit(distance, distance_idx_sorted)
  x_surrogates = ssm.randomize(x, n_rep=n_surrogate_datasets)
  surrogate_corrs = spearmanr(y, x_surrogates).flatten()
  r_stat = ss.spearmanr(x, y)
  p = nonparp(r_stat[0], surrogate_corrs)
  return [r_stat, p]


def spin_pearson(x,y):
  # x, y should be one array 
  n_surrogate_datasets = 1000

  # Note: number samples must be greater than number neighbors
  num_samples = 25
  num_neighbors = 30
  distance= dijkstra(np.array(pd.read_csv('../data/LeftParcelGeodesicDistmat.txt',
                                          header=None, delimiter=' ')), directed=False)
  distance_idx_sorted = np.argsort(distance, axis=1)
  ssm = SampledSurrogateMaps(ns=num_samples, knn=num_samples, random_state=0,resample=True)
  ssm.fit(distance, distance_idx_sorted)
  x_surrogates = ssm.randomize(x, n_rep=n_surrogate_datasets)
  surrogate_corrs = pearsonr(y, x_surrogates).flatten()
  r_stat = ss.pearsonr(x, y)
  p = nonparp(r_stat[0], surrogate_corrs)
  return [r_stat, p]

class spin_pearson_pair:
  def __init__(self,x,y1,y2):
    n_surrogate_datasets = 1000

  # Note: number samples must be greater than number neighbors
    num_samples = 40
    num_neighbors = 20
    distance = dijkstra(np.array(pd.read_csv('../data/LeftParcelGeodesicDistmat.txt',
                        header=None, delimiter=' ')), directed=False)
    distance_idx_sorted = np.argsort(distance, axis=1)
    ssm = SampledSurrogateMaps(ns=num_samples, knn=num_samples, random_state=0,resample=True)
    ssm.fit(distance, distance_idx_sorted)
    x_surrogates = ssm.randomize(x, n_rep=n_surrogate_datasets)
    self.r1_permut = pearsonr(y1, x_surrogates).flatten()
    self.r2_permut = pearsonr(y2, x_surrogates).flatten()
    r_stat1 = ss.pearsonr(x, y1)[0]
    r_stat2 = ss.pearsonr(x, y2)[0]
    self.spin_p1 = nonparp(r_stat1, self.r1_permut)
    self.spin_p2 = nonparp(r_stat2, self.r2_permut)

rng = np.random.default_rng()
def my_statistic(x, y):
  return ss.pearsonr(x, y)[0]
def res(x,y):
  a = ss.bootstrap((x, y), my_statistic, vectorized=False, paired=True,
                  random_state=rng)
  return a

def se(data):
    error = np.std(data, ddof=1) / np.sqrt(np.size(data))
    return error

def cohen_d(d1,d2):
	# calculate the size of samples
  n1, n2 = len(d1), len(d2)
	# calculate the variance of the samples
  s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
	# calculate the pooled standard deviation
  s = math.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
	# calculate the means of the samples
  u1, u2 = np.mean(d1), np.mean(d2)
	# calculate the effect size
  return (u1 - u2) / s

def permut(actual, cal, m): # x,y,n_permutations
  test_stat = ss.pearsonr(actual, cal)[0]
  perm = np.array([np.random.permutation(actual) for _ in range(m)])
  naive_corrs = [ss.pearsonr(perm[i], cal)[0] for i in range(m)]
  p = np.sum(np.abs(naive_corrs) > abs(test_stat)) / m
  return p

def partial_r(A, B, C): # A, B, C should be the 1-d array
    AB = ss.pearsonr(A,B)[0]
    AC = ss.pearsonr(A,C)[0]
    BC = ss.pearsonr(B,C)[0]
    partial = AB - AC*BC
    deno = (1-AC*AC) * (1-BC*BC)
    r = partial/pow(deno, 0.5)
    return r

def fdr(p_vals):
    ranked_p_values = ss.rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    return fdr

def pearson_na(x,y):
  import numpy.ma as ma
  a=ma.masked_invalid(x)
  b=ma.masked_invalid(y)
  msk = (~a.mask & ~b.mask)
  rp = ss.pearsonr(a[msk],b[msk])
  return np.append(rp,len(msk))

def spearman_na(x,y):
  import numpy.ma as ma
  a=ma.masked_invalid(x)
  b=ma.masked_invalid(y)
  msk = (~a.mask & ~b.mask)
  rp = ss.spearmanr(a[msk],b[msk])
  return np.append(rp,len(msk))
    
