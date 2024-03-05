import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from scipy.signal import savgol_filter
from sklearn.manifold import TSNE


### Apply savgol filter for a wavelength filter, 
def filter_segment(features_noWtab, order=1,der= False):
    #features_noWtab: Segment of the signal
    #order: Order of the savgol filter
    #der: If with first derivative
    
#     part1 = features_noWtab.loc[:,indx]
    part1 = features_noWtab.copy()
    if (der):
        fr1 = savgol_filter(part1, 65, 1,deriv=1)
    else:
        fr1 = savgol_filter(part1, 65, order)
    fr1 = pd.DataFrame(data=fr1, columns=part1.columns)
    return fr1


def feature_preparation(features, inval = [1351,1431, 1801, 2051], frmax=2451, order=1,der= False):
    # features: The original reflectance signal
    #order: Order of the savgol filter
    #der: If with first derivative
    
    features.columns = features.columns.astype('int')
    features[features<0] = 0   
    
    #####Substitute high values with the mean of neighbour values
    other = features.copy()
    other[other>1] = np.nan
    other = (other.fillna(method='ffill') + other.fillna(method='bfill'))/2
    other=other.interpolate(method='linear', axis=1).ffill().bfill()
    
    wt_ab = [i for i in range(inval[0],inval[1])]+[i for i in range(inval[2],inval[3])]+[i for i in range(2451,2501)] 

    features_Wtab = other.loc[:,wt_ab]
    features_noWtab=other.drop(wt_ab,axis=1)
    
    fr1 = filter_segment(features_noWtab.loc[:,:inval[0]-1], order = order, der = der)
    fr2 = filter_segment(features_noWtab.loc[:,inval[1]:inval[2]-1], order = order,der = der)
    fr3 = filter_segment(features_noWtab.loc[:,inval[3]:frmax], order = order,der = der)    
    
    
    inter = pd.concat([fr1,fr2,fr3], axis=1, join='inner')
    inter[inter<0]=0
    
    return inter

## Plot spectra signal for all saomles 
def plot_fig(features, save=False, file=None, figsize=(15, 10)):
    # features: The original reflectance signal
    
    plt.figure(figsize=figsize)
    plt.plot(features.T)
    plt.ylim(0, features.max().max())
    if (save):
        plt.savefig(file + '.pdf', bbox_inches = 'tight', dpi = 1000)
        plt.savefig(file + '.svg', bbox_inches = 'tight', dpi = 1000)
    plt.show()


## Plot the statsitical distribution of reflectance values across the full range##
def spectraRange(inter, inval = [1351,1431, 1801, 2051], save=False, file=None, figsize=(10, 5)):
    # inter: Reflectance signal    
    import scipy.stats as st

    inter.columns = inter.columns.astype('int')
    inter = inter.loc[:,400:].copy()
    
    wt_ab = [i for i in range(inval[0],inval[1])]+[i for i in range(inval[2],inval[3])]+[i for i in range(2451,2501)] 

    for j in wt_ab:
        inter[j] = np.nan

    inter = inter.T.sort_index().T
    
    smooth_path = inter.describe().T['mean']#np.mean(inter)          #sample mean
    min_s = inter.quantile(0.01)
    max_s = inter.quantile(0.99)
    
    plt.rc('font', size=15)
    plt.figure(figsize=figsize)

    # Plotting:
    plt.plot(min_s, linewidth=1, label='1-99% quantile',linestyle='--',c='black')
    plt.plot(max_s, linewidth=1,linestyle='--',c='black')
    plt.plot(smooth_path.T, linewidth=2, label='Mean',c='black')  # mean curve.

    plt.ylabel('Reflectance')
    plt.xlabel('Wavelength (nm)')

    plt.legend()

    if (save):
        plt.savefig('{}.svg'.format(file), bbox_inches='tight', dpi=300)
        plt.savefig('{}.png'.format(file), bbox_inches='tight', dpi=300)


def meanLandCoverRange(inter, w, inval = [1351,1431, 1801, 2051], save=False, file=None):
    # inter: Reflectance signal
    # w: Metadata
    
    plt.rcParams['lines.markersize'] = 0.5

    plt.rc('font', size=11)
    fig , axes= plt.subplots(nrows=2, ncols=1 ,dpi=500,constrained_layout=True)

    lc = ['Tundra', 'Forest', 'Crops', 'Grassland', 'Shrubland','Mix']
    colors = {'Tundra': '#A64B00', 
    'Forest': '#008500', 
    'Crops': '#85004B', 
    'Grassland': 'lightgreen', 
    'Shrubland': '#A60000',
    'Mix': '#FF0000'}

    # wt_ab = [i for i in range(1351, 1501)] + [i for i in range(1800, 1951)] + [i for i in range(2450, 2501)]
    wt_ab = [i for i in range(inval[0],inval[1])]+[i for i in range(inval[2],inval[3])]+[i for i in range(2451,2501)] 
    inter.columns = inter.columns.astype('int')
    groups = pd.concat([w,inter],axis=1).groupby(['LandCover'])

    for i in lc:
        group = groups.get_group(i).loc[:,400:]        
        
        for j in wt_ab:
            group[j] = np.nan

        des = group.describe().T.sort_index()
        var = (des['std']/des['mean'])
        des['mean'].sort_index().plot(color=colors[i],ax=axes[0], title = '(a)')
        var.plot(color=colors[i],ax=axes[1], title = '(b)')

    plt.legend(lc,bbox_to_anchor=(0.8,2.8),loc='upper left')
    
    if (save):
        plt.savefig('{}.svg'.format(file), bbox_inches='tight', dpi=500)
        plt.savefig('{}.pdf'.format(file), bbox_inches='tight', dpi=500)

## Plot the correlation between trait values ##
def plot_corr(a, method,save=False, file=None, b = None, annot=False, figsize=(20, 15)):
    plt.rc('font', size=8)
    
    if(b is not None):
        con = pd.concat([a, b], axis=1, join='inner')
        C_mat = con.corr()#.fillna(0)
    else:
        C_mat = a.corr(method = method)#.fillna(0)
    
    # Configure a custom diverging colormap #
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    mask = np.triu(np.ones_like(C_mat, dtype=bool))

    with sns.axes_style("white"):
        f, ax = plt.subplots(figsize= figsize)
        ax = sns.heatmap(C_mat, mask=mask, cmap=cmap,annot=annot, vmin=-0.8, vmax=0.8,square=True)
    if(save):
        f.savefig(file+ '.pdf', bbox_inches='tight', dpi=300)
        f.savefig(file+ '.svg', bbox_inches='tight', dpi=300)

## Color palette ##
def get_cmap(n, name='hsv'):
    return plt.cm.get_cmap(name, n)

## Plot projection of the reflectance signal of the complete data set
def pca_plot(inter, w , sp = ['LandCover','dataset'], save=False, file=None, tsne=False,figsize=(20, 15)):

    plt.rcParams['lines.markersize'] = 10
    plt.rc('font', size=15)
    fig, ax = plt.subplots(figsize=figsize)
    
    if (tsne):
        X_embedded = TSNE(n_components=2).fit_transform(inter.iloc[:, :].values)
        x = X_embedded[:, 0]
        y = X_embedded[:, 1]
    else:
        pca = PCA(n_components=2)
        dt = pca.fit_transform(inter.values)
        x = dt[:, 0]
        y = dt[:, 1]

    markers = {'Tundra': ".", 'Forest': ",", 'Crops': "o", 'Grassland': "v", 'Shrubland': "^", 'Mix': "<"}
    merge = pd.concat([w, pd.DataFrame([x,y], index=['x','y']).T], axis=1)

    for a, gp in merge.groupby(sp):
        if (sp == 'LandCover'):
            colors = {'Tundra': '#A64B00', 
    'Forest': '#008500', 
    'Crops': '#85004B', 
    'Grassland': 'lightgreen', 
    'Shrubland': '#A60000',
    'Mix': '#FF0000'}
    #         ax.scatter(gp.x, gp.y, marker = markers[gp[sp].unique()[0]], label = gp[sp].unique()[0], color = colors[gp[sp].unique()[0]])
            ax.scatter(gp.x, gp.y, label = gp[sp].unique()[0], color = colors[gp[sp].unique()[0]])
        elif(sp == 'dataset'):
            ax.scatter(gp.x, gp.y, label = 'DB'+str(gp[sp].unique()[0]))        
        else:
            ax.scatter(gp.x, gp.y,  marker = markers[gp['LandCover'].unique()[0]], label = 'DB'+str(gp['dataset'].unique()[0]) +': ' +gp['LandCover'].unique()[0])

    ax.legend()
    if (save):
        plt.savefig(file+ '.pdf', bbox_inches='tight', dpi=1000)
        plt.savefig(file+ '.svg', bbox_inches='tight', dpi=1000)

    plt.show()

## Bar plot with the number of non-null samples per Trait ##
def labels_bar(labels, Traits, figsize=(15, 5)):
    ax = labels.describe().T['count'].reset_index(name='count').sort_values(['count'], ascending=False).plot(kind='bar',
                                                                                                             figsize=figsize)
    ax.bar_label(ax.containers[0])

    plt.show()