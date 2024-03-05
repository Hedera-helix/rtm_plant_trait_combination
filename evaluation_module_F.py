from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error,mean_squared_error, r2_score 
from sklearn.preprocessing import RobustScaler,MinMaxScaler,PowerTransformer,StandardScaler
from sklearn.pipeline import Pipeline

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import math
from pickle import dump,load

from model_module_F import samp_w
from feature_module_F import plot_fig
from scipy import stats

import sys
import os

def collectPredictions(paths, name=''):
#     data = pd.read_csv(paths[0] + '{}_{}.csv'.format(name, 0)).drop(['Unnamed: 0'],axis=1)
    data = pd.DataFrame()

    for dir_n in paths:
        sc = pd.read_csv(dir_n + '{}.csv'.format(name)).drop(['Unnamed: 0'],axis=1)
        data = pd.concat([data,sc],axis=0)

    return data.reset_index(drop=True)

def all_scores(test_tr,Traits,obs_pf, pred_df,samp_w_ts=None, method = None, save =False, dir_n = None):
    r2_tab = []
    RMSE_tab = []
    nrmse_tab = []
    mae_tab = []
    b_tab = []

    for j in test_tr:

        f = pred_df[j+ ' Predictions'].reset_index(drop=True) # + ' Predictions'
        y = obs_pf[j].reset_index(drop=True)

        idx = np.union1d(f[f.isna()].index,y[y.isna()].index)

        f.drop(idx, axis = 0, inplace=True)
        y.drop(idx, axis = 0, inplace=True)
        
        
        if (y.notnull().sum()):
            if (samp_w_ts is not None):
                we = pd.DataFrame(samp_w_ts).loc[f.index,:]
            else:
                we = None

            if (we is not None) and (we.sum().sum() !=0):
                r2_tab.append(r2_score(y,f,sample_weight= we))

                RMSE=math.sqrt(mean_squared_error(y,f,sample_weight= we))
                RMSE_tab.append(RMSE)
                nrmse_tab.append((RMSE*100)/(np.nanquantile(np.array(y),0.99) - np.nanquantile(np.array(y),0.01)))

                mae_tab.append(mean_absolute_error(y,f,sample_weight= we))

                bias=np.sum(np.array(y)-np.array(f))/len(f)
                b_tab.append(bias)
            else:
                r2_tab.append(r2_score(y,f))

                RMSE=math.sqrt(mean_squared_error(y,f))
                RMSE_tab.append(RMSE)
                nrmse_tab.append((RMSE*100)/(np.nanquantile(np.array(y),0.99) - np.nanquantile(np.array(y),0.01)))

                mae_tab.append(mean_absolute_error(y,f))

                bias=np.sum(np.array(y)-np.array(f))/len(f)
                b_tab.append(bias)
        else:
            r2_tab.append(np.nan)
            RMSE_tab.append(np.nan)
            nrmse_tab.append(np.nan)
            mae_tab.append(np.nan)
            b_tab.append(np.nan)
            pass        

    test = pd.DataFrame([r2_tab, RMSE_tab, nrmse_tab,mae_tab,b_tab], columns= test_tr[:len(test_tr)], index=['r2_score','RMSE','nRMSE (%)','MAE','Bias'])
    if(save):
        test.to_csv(dir_n + 'CV_metrics_{}.csv'.format(method))
    return test



def scatterPlot(obs, preds, meta, Traits, test_tr, test, sp = None, method = None, save =False, dir_n = None, figsize=(7.48, 15)):

    if sp is not None:
        n = len(meta[sp].unique())
        cmap = sns.color_palette('colorblind', n)

    plt.rc('font', size=7)
    plt.rcParams['lines.markersize'] = 3
    plt.rcParams['lines.linewidth'] = 0.5

    raws = round((len(test_tr)+1)/2)
    count=1
    ax1 = plt.subplots(figsize= figsize, dpi=300 ,sharex=True, sharey=True,constrained_layout=True)

    for j in range(len(test_tr)):
        f = preds.loc[:,Traits[Traits.index(test_tr[j])]+' Predictions']
        # y = obs.loc[f.index,Traits[Traits.index(test_tr[j])]]
        y = obs.loc[:,Traits[Traits.index(test_tr[j])]]

        idx = np.union1d(f[f.isna()].index,y[y.isna()].index)

        f.drop(idx, axis = 0, inplace=True)
        y.drop(idx, axis = 0, inplace=True)

        m = meta.drop(idx, axis = 0)

        plt.subplot(raws,4,count)

        lim = max(f.max(),y.max())  #min(f.quantile(0.99),y.quantile(0.99)) min(f.max(),y.max()) max(f.quantile(0.99),y.quantile(0.99))
        ax1 = sns.lineplot(x=(0,lim), y=(0,lim), color='black',legend='full', linestyle='dashed')          
        # ax1.set_xlim(0, lim)
        # ax1.set_ylim(ax1.get_xlim())
        
#         label_format = '{:,.3f}'
#         ticks_loc = ax1.get_xticks().tolist()
#         ax1.set_yticks(ax1.get_xticks().tolist())
#         ax1.set_yticklabels([label_format.format(x) for x in ticks_loc])
        
#         ax1.set_xticks(ax1.get_xticks().tolist())
#         ax1.set_xticklabels([label_format.format(x) for x in ticks_loc])
        
        # ax1.set_aspect('equal', 'datalim')
        ax1.set_aspect(1.0/ax1.get_data_ratio(), adjustable='box')
        # ax1.set_box_aspect(1)

        y.name = y.name +' Observations'

        slope, intercept, r_value, p_value, std_err = stats.linregress(f,y)

        ax1 = sns.regplot(x= f, y=y, color='b', fit_reg= True, scatter_kws={"color": "beige"})

        if sp is not None:
            groups = pd.concat([f,y,m], axis=1).groupby(sp)

            for name, group in groups:
                sns.regplot(x = group[test_tr[j]+' Predictions'], y = group[test_tr[j]+' Observations'] ,fit_reg= False, ci=False, label= sp + ' {}'.format(group[sp].unique()[0]),ax = ax1, scatter_kws={ "color": cmap[list(meta[sp].unique()).index(name)], 'alpha':0.3})
        
        else:
            sns.regplot(x = f, y = y ,fit_reg= False, ci=False, scatter_kws={"color": "blue", 'alpha':0.3})

        if j==0:
            handles, labels = ax1.get_legend_handles_labels()
            
        # ax1.set_xlim(0, lim)
        # ax1.set_ylim(ax1.get_xlim())
        # ax1.set_aspect(1.0/ax1.get_data_ratio(), adjustable='box')

        ann = 'y = {0:.2f}x+{1:.2f} \n R² = {2:.2f} \n nRMSE = {3:.2f}'.format(slope,intercept,test.loc['r2_score',test_tr][j],test.loc['nRMSE (%)',test_tr][j])

        ax1.annotate(ann,
            xy=(0.5,0.01),
            xycoords='axes fraction',
            horizontalalignment='left',
            verticalalignment='bottom',size=5.75)

        ann = test_tr[j]
        ax1.set_title(ann, y=1.1, pad=-5, fontdict = {'fontsize':5,
     'horizontalalignment': 'center', 'fontweight':'bold'})

        plt.xlabel(" ")
        plt.ylabel(" ")

        count+=1

    plt.legend(handles, labels, loc='best')

    if(save):
        plt.savefig(dir_n + "Traits_scatter_plot_allDataset_{}.pdf".format(method),bbox_inches = 'tight', dpi = 300)
        plt.savefig(dir_n + "Traits_scatter_plot_allDataset_{}.svg".format(method),bbox_inches = 'tight', dpi = 300)

# def scatterPlot(obs, preds, meta, Traits, test_tr, test, method = None, save =False, dir_n = None,figsize=(7.48, 15)):
#     plt.rc('font', size=7)
#     plt.rcParams['lines.markersize'] = 1
#     plt.rcParams['lines.linewidth'] = 0.5

#     raws = round((len(test_tr)+1)/2)
#     count=1
#     ax = plt.subplots(figsize=figsize, dpi=1000 ,sharex=True,constrained_layout=True)

#     for j in range(len(test_tr)):
#         f = preds.loc[:,Traits[Traits.index(test_tr[j])]+' Predictions'] #+' Predictions'
#         y = obs.loc[:,Traits[Traits.index(test_tr[j])]]

#         idx = np.union1d(f[f.isna()].index,y[y.isna()].index)

#         f.drop(idx, axis = 0, inplace=True)
#         y.drop(idx, axis = 0, inplace=True)

#         # we = pd.DataFrame(weights).drop(idx, axis = 0)
#         m = meta.drop(idx, axis = 0)

#         plt.subplot(raws,4,count)


#         lims = [0, max(max(np.array(f).reshape(1,-1)[0]),max(np.array(y).reshape(1,-1)[0]))] 

#         plt.xlim(lims)
#         plt.ylim(lims)    
#         _ = plt.plot(lims, lims,label='Identity line')    

#         y.name = y.name +' Observations'
#     #     f.name = f.name +' Predictions'

#         slope, intercept, r_value, p_value, std_err = stats.linregress(f,y)

#         sns.regplot(x= f, color='b', y=y,fit_reg= True,scatter_kws={"color": "beige"})
#         groups = pd.concat([f,y,m], axis=1).groupby("LandCover")
#         colors = {'Tundra': '#A64B00', 
#         'Forest': '#008500', 
#         'Crops': '#85004B', 
#         'Grassland': 'lightgreen', 
#         'Shrubland': '#A60000',
#         'Mix': '#FF0000'}

#         for name, group in groups:
#             plt.scatter(group[test_tr[j]+' Predictions'], group[test_tr[j]+' Observations'], c=colors[group['LandCover'].unique()[0]], alpha = 0.3, label= group['LandCover'].unique()[0])

#         plt.xlabel(" ")
#         plt.ylabel(" ")
#         ann = 'y = {0:.2f}x+{1:.2f} \n R² = {2:.2f} \n nRMSE = {3:.2f}'.format(slope,intercept,test.loc['r2_score',test_tr][j],test.loc['nRMSE (%)',test_tr][j])
#         plt.annotate(ann,
#                 xy=(max(max(np.array(f).reshape(1,-1)[0]),max(np.array(y).reshape(1,-1)[0]))/2.25,min(max(np.array(f).reshape(1,-1)[0]),max(np.array(y).reshape(1,-1)[0])/30)),
#                 xycoords='data',
#                 horizontalalignment='left',
#                 verticalalignment='bottom',size=5.75)
#         ann = test_tr[j]
#         plt.annotate(ann,
#                 xy=(0,max(max(np.array(f).reshape(1,-1)[0]),max(np.array(y).reshape(1,-1)[0]))),
#                 xycoords='data',
#                 horizontalalignment='left',
#                 verticalalignment='top',size=4.75, fontweight='bold')

#         plt.axis('square')
#         # plt.legend(fontsize=7)

#         count+=1
        
#     if(save):
#         plt.savefig(dir_n + "Traits_scatter_plot_allDataset_{}.pdf".format(method),bbox_inches = 'tight', dpi = 300)
#         plt.savefig(dir_n + "Traits_scatter_plot_allDataset_{}.svg".format(method),bbox_inches = 'tight', dpi = 300)

        
def histPlot(obs, preds, meta, Traits, test_tr, method = None, save =False, dir_n = None,figsize=(15,40)):
    
    plt.rc('font', size=10)
    plt.rcParams['lines.markersize'] = 1
    plt.rcParams['lines.linewidth'] = 0.5
    
    a = round((len(test_tr)+1)/2)  # number of rows
    b = 2  # number of columns
    c = 1  # initialize plot counter

    fig = plt.figure(figsize=figsize)

    for j in range(len(test_tr)):

        f = preds.loc[:,Traits[Traits.index(test_tr[j])]+' Predictions']
        y = obs.loc[:,Traits[Traits.index(test_tr[j])]]

        idx = np.union1d(f[f.isna()].index,y[y.isna()].index)

        f.drop(idx, axis = 0, inplace=True)
        y.drop(idx, axis = 0, inplace=True)

        # we = pd.DataFrame(weights).drop(idx, axis = 0)
        m = meta.drop(idx, axis = 0)    

        test = pd.DataFrame(np.array([y,f]).T , columns=['Observations','Predictions']) 

        plt.subplot(a, b, c)
        plt.xlabel(test_tr[j])
        sns.histplot(test, alpha=0.3 ,element="step", label = test_tr[j],stat="probability")

        c = c + 1
        
    if(save):
        plt.savefig(dir_n + "Traits_hist_plot_allDataset_{}.pdf".format(method),bbox_inches = 'tight', dpi = 600)
        plt.savefig(dir_n + "Traits_hist_plot_allDataset_{}.svg".format(method),bbox_inches = 'tight', dpi = 600)

def visualize_duplicates(gap_fil, Traits, tr, train_x, save=False,dir_n=None):
    #### Plot the maximum duplicated values of trait tr###
    df = gap_fil[gap_fil[Traits[tr]].duplicated(keep=False)][Traits]
    if (df[Traits[tr]].sum()!=0):
        df1 = (df.groupby(Traits[tr])
           .apply(lambda x: tuple(x.index))       
           .reset_index(name='idx'))    
        t = np.argmax([len(df1['idx'][i]) for i in df1.index])
        gap_fil.loc[df1.loc[t, 'idx'],train_x.columns].T.plot(figsize=(15, 10))
        plt.title('Duplicated values ({}) of {}'.format(df1.loc[t, Traits[tr]], Traits[tr]))
        if (save):
            plt.savefig(dir_n + '/plot_Duplicatedvalues_Tr{}.pdf'.format(tr), bbox_inches = 'tight', dpi = 1000)
            plt.savefig(dir_n + '/plot_Duplicatedvalues_Tr{}.svg'.format(tr), bbox_inches = 'tight', dpi = 1000)
    else:
        pass

def dropped_spectra(db_train, gap_fil, save=False, dir_n=None):
    db = db_train.loc[:,'400':]
    db.columns = db.columns.astype('int')
    drop = db.loc[np.array(db_train.index.drop(gap_fil.index)),:]
    if(save):
        gap_fil.index.to_csv(dir_n + '/droppedSpectra.csv')
        plot_fig(drop, save=True, file= dir_n + '/plot_droppedSpectra')
    
def cleanana(score):
    r = score
    r[r<0] = 0
    r[r>1] = 5
    return r