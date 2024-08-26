import numpy as np
import hcp_utils as hcp
import seaborn as sns
import matplotlib.pyplot as plt

network_name = np.array(['Vis1', 'Vis2', 'SMN', 'CON', 'DAN', 'Lan', 
                         'FPN', 'Aud', 'DMN', 'PMN', 'VMN', 'OAN'])
                         
rgb = np.array(list(hcp.ca_network['rgba'].values())[1:])

def plot_t_single(t_g123, thres, filename):
  rank = t_g123.argsort()
  sns.set_context("paper", font_scale = 1.5)
  fig, ax = plt.subplots(figsize=(6,3))
  ax.bar(network_name[rank], t_g123[rank], color=rgb[rank])
  if t_g123.min() < 0:
    ax.set_ylim([-thres,thres])
  else:
    ax.set_ylim([0,thres])
  ax.set_yticks([])
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  ax.spines['left'].set_visible(False)
  fig.tight_layout()
  fig.savefig(filename, dpi=300, transparent=True)