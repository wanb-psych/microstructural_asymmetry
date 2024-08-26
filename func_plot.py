import numpy as np
import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
import seaborn as sns
from brainspace.plotting import plot_hemispheres
from brainspace.utils.parcellation import map_to_labels
import hcp_utils as hcp
from brainsmash.mapgen.stats import spearmanr
from brainsmash.mapgen.stats import pearsonr
from brainsmash.mapgen.stats import nonparp
from brainspace.null_models import SampledSurrogateMaps
from scipy.sparse.csgraph import dijkstra
from brainspace.datasets import load_conte69
from scipy.interpolate import make_interp_spline
from scipy.io import loadmat
from matplotlib.colors import ListedColormap
import math
# lh, rh = load_conte69()

c_color = np.vstack((list(hcp.ca_network['rgba'].values())))[1:]
cmap_ca = ListedColormap(c_color)
cmap_ca_black = ListedColormap(np.concatenate(([[0,0,0,1]],c_color)))
cmap_ca_white = ListedColormap(np.concatenate(([[1,1,1,1]],c_color)))

ct_color = [[1,1,1],
            [0.9961, 0.9961, 0.9961],
            [0.9844, 0.8945, 0.9844],
            [0.9883, 0.8242, 0.7812],
            [0.8516, 0.7734, 0.5977],
            [0.6680, 0.7266, 0.6328],
            [0.5430, 0.6523, 0.6875],
            [0.4961, 0.5469, 0.6719]]
cmap_ct = ListedColormap(ct_color)


def plot_surface_ll(lh, rh,
                    data, 
                    size,
                    cmap,
                    color_range,
                    filename):
  plot_hemispheres(lh, rh, array_name = data, nan_color = (0,0,0,1),size = size,
                   cmap = cmap, color_bar = True, color_range=color_range,
                   interactive = False, zoom = 1.5, embed_nb = True, transparent_bg=True,
                   cb__labelTextProperty={"fontSize": 24}, screenshot=True, filename=filename)
  fig = plot_hemispheres(lh, rh, array_name = data, nan_color = (0,0,0,1),size = size,
                         cmap = cmap, color_bar = True, color_range=color_range,
                         cb__labelTextProperty={"fontSize": 24}, interactive = False, zoom = 1.5, embed_nb = True)
  return fig

def plot_surface_lr(lh, rh,
                    data, 
                    size,
                    cmap,
                    color_range,
                    filename):
  plot_hemispheres(lh, rh, array_name = data, nan_color = (0,0,0,1),size = size,
                   cmap = cmap, color_bar = True, color_range=color_range,
                   interactive = False, zoom = 1.5, embed_nb = True, transparent_bg=True,
                   cb__labelTextProperty={"fontSize": 24}, screenshot=True, filename=filename)
  fig = plot_hemispheres(lh, rh, array_name = data, nan_color = (0,0,0,1),size = size,
                         cmap = cmap, color_bar = True, color_range=color_range,
                         cb__labelTextProperty={"fontSize": 24}, interactive = False, zoom = 1.5, embed_nb = True)
  return fig

network_name = np.array(['Vis1', 'Vis2', 'SMN', 'CON', 'DAN', 'LAN', 
                         'FPN', 'AUD', 'DMN', 'PMN', 'VMN', 'OAN'])
                         
rgb = np.array(list(hcp.ca_network['rgba'].values())[1:])

def plot_t_single(t_g123, thres, bar, ax, x_rota=0, x_name=network_name, color=rgb,
                  ascend=True, y_axis=False):
  rank = t_g123.argsort()
  if ascend==False:
    rank = rank[::-1] 
  if len(np.array(color).shape) > 1:
    color = color[rank]
  if type(bar) == type(None):
    ax.bar(x_name[rank], t_g123[rank], color=color)
  else:
    ax.bar(x_name[rank], t_g123[rank], yerr= bar[rank], color=color, align='center',
           error_kw=dict(ecolor='gray', zorder=0, lw=2, capsize=5, capthick=2))
  ax.set_xticklabels(network_name[rank], rotation = x_rota)
  if t_g123.min() < 0:
    ax.set_ylim([-thres,thres])
    ax.set_yticks([-thres,0,thres])
  else:
    ax.set_ylim([0,thres])
    ax.set_yticks([0,thres])
  if y_axis==False:
    ax.set_yticks([])
    ax.spines['right'].set_visible(False)
  
  ax.spines['top'].set_visible(False)
  ax.spines['left'].set_visible(False)

def plot_t_single_h(t_g123, thres, bar, ax,y_name=network_name, color=rgb,
                    ascend=True, x_axis=False):
  rank = t_g123.argsort()
  if len(np.array(color).shape) > 1:
    color = color[rank]
  if ascend==False:
    rank = rank[::-1]
  if type(bar) == type(None):
      ax.barh(y=y_name[rank], width=t_g123[rank], color=color)
  else:
      ax.barh(y=y_name[rank], width=t_g123[rank], xerr= bar[rank], color=color, align='center',
              error_kw=dict(ecolor='gray', zorder=0, lw=2, capsize=5, capthick=2))

  if t_g123.min() < 0:
    ax.set_xlim([-thres,thres])
    ax.set_xticks([-thres,0,thres])
  else:
    ax.set_xlim([0,thres])
    ax.set_xticks([0,thres])
  if x_axis==False:
    ax.set_xticks([])
    ax.spines['bottom'].set_visible(False)
  
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)

def ml_r_displot(data,legend='off'):
  sns.set_context("paper", font_scale = 2.5)
  fig, ax = plt.subplots(1, figsize=(8,6))
  L1 = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0']
  for i in range(10):
    sns.distplot(data[i], hist = False, label='L1_ratio='+L1[i], ax=ax)
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  ax.set_xlabel('Pearson $\it{r}$')
  if legend=='on':
    ax.legend(loc='upper left', markerscale=0.7, scatterpoints=1, fontsize=14, frameon=False)
  fig.tight_layout()


from matplotlib.colors import LinearSegmentedColormap, ListedColormap

#__all__ = ['parula', 'justine', 'dinosaur']

parula = LinearSegmentedColormap.from_list(
    'parula',
    ['#352A87', '#363093', '#3637A0', '#353DAD', '#3243BA', '#2C4AC7',
     '#2053D4', '#0F5CDD', '#0363E1', '#0268E1', '#046DE0', '#0871DE',
     '#0D75DC', '#1079DA', '#127DD8', '#1481D6', '#1485D4', '#1389D3',
     '#108ED2', '#0C93D2', '#0998D1', '#079CCF', '#06A0CD', '#06A4CA',
     '#06A7C6', '#07A9C2', '#0AACBE', '#0FAEB9', '#15B1B4', '#1DB3AF',
     '#25B5A9', '#2EB7A4', '#38B99E', '#42BB98', '#4DBC92', '#59BD8C',
     '#65BE86', '#71BF80', '#7CBF7B', '#87BF77', '#92BF73', '#9CBF6F',
     '#A5BE6B', '#AEBE67', '#B7BD64', '#C0BC60', '#C8BC5D', '#D1BB59',
     '#D9BA56', '#E1B952', '#E9B94E', '#F1B94A', '#F8BB44', '#FDBE3D',
     '#FFC337', '#FEC832', '#FCCE2E', '#FAD32A', '#F7D826', '#F5DE21',
     '#F5E41D', '#F5EB18', '#F6F313', '#F9FB0E']
)

justine = ListedColormap(
    ['#3855A5', '#3857A6', '#3958A8', '#395AA9', '#395BAB', '#3A5DAC',
     '#3A5EAD', '#3A60AF', '#3B61B0', '#3B63B1', '#3B65B3', '#3C66B4',
     '#3C68B6', '#3D69B7', '#3D6BB8', '#3D6CBA', '#3E6EBB', '#3E6FBC',
     '#3E71BE', '#3F73BF', '#3F74C1', '#3F76C2', '#4077C3', '#4079C5',
     '#407AC6', '#417CC8', '#417DC9', '#417FCA', '#4281CC', '#4282CD',
     '#4284CE', '#4385D0', '#4387D1', '#4488D3', '#448AD4', '#448BD5',
     '#458DD7', '#458FD8', '#4590D9', '#4692DB', '#4693DC', '#4695DE',
     '#4796DF', '#4798E0', '#4799E2', '#489BE3', '#489DE5', '#489EE6',
     '#49A0E7', '#49A1E9', '#49A3EA', '#4AA4EB', '#4AA6ED', '#4BA7EE',
     '#4BA9F0', '#4BABF1', '#4CACF2', '#4CAEF4', '#4CAFF5', '#4DB1F6',
     '#4DB2F8', '#4DB4F9', '#4EB5FB', '#4EB7FC', '#51B8FC', '#54B9FC',
     '#56BAFC', '#59BBFC', '#5CBDFC', '#5FBEFC', '#61BFFC', '#64C0FC',
     '#67C1FC', '#6AC2FC', '#6CC3FD', '#6FC5FD', '#72C6FD', '#75C7FD',
     '#77C8FD', '#7AC9FD', '#7DCAFD', '#80CBFD', '#83CCFD', '#85CDFD',
     '#88CFFD', '#8BD0FD', '#8ED1FD', '#90D2FD', '#93D3FD', '#96D4FD',
     '#99D5FD', '#9BD7FD', '#9ED8FD', '#A1D9FD', '#A4DAFD', '#A6DBFE',
     '#A9DCFE', '#ACDDFE', '#AFDEFE', '#B2DFFE', '#B4E1FE', '#B7E2FE',
     '#BAE3FE', '#BDE4FE', '#BFE5FE', '#C2E6FE', '#C5E7FE', '#C8E8FE',
     '#CAEAFE', '#CDEBFE', '#D0ECFE', '#D3EDFE', '#D6EEFE', '#D8EFFE',
     '#DBF0FE', '#DEF2FE', '#E1F3FE', '#E3F4FF', '#E6F5FF', '#E9F6FF',
     '#ECF7FF', '#EEF8FF', '#F1F9FF', '#F4FAFF', '#F7FCFF', '#F9FDFF',
     '#FCFEFF', '#FFFFFF', '#FFFEFC', '#FFFDF9', '#FFFCF7', '#FFFBF4',
     '#FFF9F1', '#FFF8EE', '#FFF7EC', '#FFF6E9', '#FFF5E6', '#FFF4E3',
     '#FFF3E1', '#FFF2DE', '#FFF0DB', '#FFEFD8', '#FFEED6', '#FFEDD3',
     '#FFECD0', '#FFEBCD', '#FFEACB', '#FFE9C8', '#FFE7C5', '#FFE6C2',
     '#FFE5C0', '#FFE4BD', '#FFE3BA', '#FFE2B7', '#FFE1B5', '#FFE0B2',
     '#FFDFAF', '#FFDDAC', '#FFDCAA', '#FFDBA7', '#FFD9A3', '#FFD79F',
     '#FFD59C', '#FFD398', '#FFD194', '#FECF90', '#FECE8C', '#FECC88',
     '#FECA85', '#FEC881', '#FEC67D', '#FEC479', '#FEC275', '#FEC072',
     '#FEBE6E', '#FEBC6A', '#FEBA66', '#FDB862', '#FDB65F', '#FDB45B',
     '#FDB257', '#FDB053', '#FDAE4F', '#FDAC4C', '#FDAA48', '#FDA844',
     '#FDA740', '#FDA53C', '#FDA339', '#FCA135', '#FC9F31', '#FC9D2D',
     '#FC9B29', '#FC9926', '#FC9722', '#FC951E', '#FB921E', '#FB901F',
     '#FA8D1F', '#F98A1F', '#F9871F', '#F88520', '#F88220', '#F77F20',
     '#F67D21', '#F67A21', '#F57721', '#F47521', '#F47222', '#F36F22',
     '#F26C22', '#F26A22', '#F16723', '#F16423', '#F06223', '#EF5F24',
     '#EF5C24', '#EE5A24', '#ED5724', '#ED5425', '#EC5125', '#EB4F25',
     '#EB4C26', '#EA4926', '#EA4726', '#E94426', '#E84127', '#E83E27',
     '#E73C27', '#E63928', '#E63628', '#E53428', '#E43128', '#E42E29',
     '#E32C29', '#E22929', '#E22629', '#DF2528', '#DD2427', '#DB2326',
     '#D82225', '#D62123', '#D42022', '#D11F21', '#CF1E20', '#CD1D1E',
     '#CA1C1D', '#C81A1C', '#C6191B', '#C31819', '#C11718', '#BE1617',
     '#BC1516', '#BA1414', '#B71313', '#B51212'],
    'justine'
)

# https://doi.org/10.1038/s41586-022-04770-6
dinosaur = LinearSegmentedColormap.from_list(
    'dinosaur',
    ['#02B2CE', '#0DB3C4', '#18B4BB', '#24B6B2', '#2FB7A9', '#3AB8A0',
     '#46BA97', '#51BB8E', '#5CBC85', '#68BE7C', '#73BF73', '#7EC06A',
     '#8AC261', '#95C358', '#A1C44F', '#ACC645', '#B7C73C', '#C3C833',
     '#CECA2A', '#D9CB21', '#E5CC18', '#F0CE0F', '#FBCF06', '#FECD04',
     '#FDC805', '#FDC405', '#FCC006', '#FBBC07', '#FBB807', '#FAB408',
     '#FAB009', '#F9AC09', '#F8A80A', '#F8A40B', '#F79F0C', '#F69B0C',
     '#F6970D', '#F5930E', '#F48F0E', '#F48B0F', '#F38710', '#F38310',
     '#F27F11', '#F17B12', '#F17612', '#F07213', '#EF6E14', '#EF6A14',
     '#EE6615', '#ED6216', '#ED5E17', '#EC5A17', '#EC5618', '#EB5219',
     '#EA4D19', '#EA491A', '#E9451B', '#E8411B', '#E83D1C', '#E7391D',
     '#E6351D', '#E6311E', '#E52D1F', '#E52920']
)

# colors I like
red=[244/256, 81/256, 30/256]
blue=[37/256, 150/256, 190/256]