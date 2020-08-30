# generate metadata and analysis of SNR in single vs cyclic and TMA replicates
# date: 2020-05-23
# author: engje
# language Python 3.6
# license: GPL>=v3

#libraries 
import os
import pandas as pd
import numpy as np
from cmIF import preprocess, process, analyze, mpimage, cmif
import shutil
#import importlib
import skimage
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

#paths
codedir = '/home/groups/graylab_share/OMERO.rdsStore/engje/Data/20191104_ValidationStudies'
# 1 get exposure time
'''
s_type = 'czi'

d_process = {
    '44290':f'{codedir}/Images/{s_type}/44290',
    #'BM-Her2N75-15':f'{codedir}/Images/{s_type}/BM-Her2N75-15',
    #'BM-Her2N75-17':f'{codedir}/Images/{s_type}/BM-Her2N75-17',
    #'BM-Her2N75-18':f'{codedir}/Images/{s_type}/BM-Her2N75-18',
    #'4165NPanc':f'{codedir}/Images/{s_type}/4165NPanc',
    #'NPPan4165':f'{codedir}/Images/{s_type}/NPPan4165',
    #'44290-146':f'{codedir}/Images/{s_type}/44290-146',
    #'B1640':f'{codedir}/Images/{s_type}/B1640',
    #'K154':f'{codedir}/Images/{s_type}/K154',
    #'K157':f'{codedir}/Images/{s_type}/K157',
    #'K175':f'{codedir}/Images/{s_type}/K175',
    #'JE-TMA-42':f'{codedir}/Images/{s_type}/JE-TMA-42',
    #'JE-TMA-43':f'{codedir}/Images/{s_type}/JE-TMA-43',
 }

for idx,(s_sample,s_path) in enumerate(d_process.items()):
    preprocess.cmif_mkdir([f'{codedir}/Metadata/{s_sample}'])

    if s_sample.find('Her2')>-1:
        os.chdir(s_path)
        df_img = mpimage.filename_dataframe()
        df_img['scene'] = [item.split('_')[-1].split('.')[0].split('-Scene-')[1] for item in df_img.index]
        df_img['slide'] = [item.split('_')[-1].split('-Scene')[0] for item in df_img.index]
        df_img['rounds'] = [item.split('_')[0] for item in df_img.index]
        cmif.exposure_times_scenes(df_img, f'{codedir}/Metadata/{s_sample}', s_path, s_end='.czi')
    elif s_sample.find('K')==0:
        df_img = cmif.parse_czi(s_path,type='r',b_scenes=True)
        cmif.exposure_times_scenes(df_img,f'{codedir}/Metadata/{s_sample}',s_path,s_end='.czi')
    elif s_sample.find('44290-146')==0:
        df_img = cmif.parse_czi(s_path,type='s',b_scenes=False)
        cmif.exposure_times(df_img,f'{codedir}/Metadata/{s_sample}',s_path)
    elif s_sample.find('JE')==0:
        df_img = cmif.parse_czi(s_path,type='s',b_scenes=True)
        cmif.exposure_times_scenes(df_img,f'{codedir}/Metadata/{s_sample}',s_path,s_end='.czi')
    else:
        df_img = cmif.parse_czi(s_path,type='r',b_scenes=False)
        cmif.exposure_times(df_img,f'{codedir}/Metadata/{s_sample}',s_path)
'''
# 2 parse tif file names (marker metadata)
'''
s_type = 'tiff'
d_process = {
    '44290':f'{codedir}/Images/{s_type}/44290',
    #'BM-Her2N75-15':f'{codedir}/Images/{s_type}/BM-Her2N75-15',
    #'BM-Her2N75-17':f'{codedir}/Images/{s_type}/BM-Her2N75-17',
    #'BM-Her2N75-18':f'{codedir}/Images/{s_type}/BM-Her2N75-18',
    #'4165NPanc':f'{codedir}/Images/{s_type}/4165NPanc',
    #'NPPan4165':f'{codedir}/Images/{s_type}/NPPan4165',
    #'44290-146':f'{codedir}/Images/{s_type}/44290-146',
    #'B1640':f'{codedir}/Images/{s_type}/B1640',
    #'K154':f'{codedir}/Images/{s_type}/K154',
    #'K157':f'{codedir}/Images/{s_type}/K157',
    #'K175':f'{codedir}/Images/{s_type}/K175',
    #'JE-TMA-42':f'{codedir}/Images/{s_type}/JE-TMA-42',
    #'JE-TMA-43':f'{codedir}/Images/{s_type}/JE-TMA-43',
 }

for idx,(s_sample,s_path) in enumerate(d_process.items()):
    preprocess.cmif_mkdir([f'{codedir}/Metadata/{s_sample}'])
    os.chdir(s_path)
    df_img = mpimage.parse_org(s_end = "_ORG.tif")
    df_img.sort_values(['rounds','scene','color'],inplace=True)
    df_img.to_csv(f'{codedir}/Metadata/{s_sample}/{s_sample}_TifImageDataframe.csv',header=True, index=True)
'''

# 3 combine the tif naming metadata and exposure time metadata
'''
d_process = {
    '44290':f'{codedir}/Metadata/44290',
    'BM-Her2N75-15':f'{codedir}/Metadata/BM-Her2N75-15',
    'BM-Her2N75-17':f'{codedir}/Metadata/BM-Her2N75-17',
    'BM-Her2N75-18':f'{codedir}/Metadata/BM-Her2N75-18',
    '4165NPanc':f'{codedir}/Metadata/4165NPanc',
    'NPPan4165':f'{codedir}/Metadata/NPPan4165',
    '44290-146':f'{codedir}/Metadata/44290-146',
    'B1640':f'{codedir}/Metadata/B1640',
    'K154':f'{codedir}/Metadata/K154',
    'K157':f'{codedir}/Metadata/K157',
    'K175':f'{codedir}/Metadata/K175',
    'JE-TMA-42':f'{codedir}/Metadata/JE-TMA-42',
    'JE-TMA-43':f'{codedir}/Metadata/JE-TMA-43',
 }

for idx,(s_sample,s_path) in enumerate(d_process.items()):
    df_img_exp = pd.DataFrame()
    df_img = pd.read_csv(f'{s_path}/{s_sample}_TifImageDataframe.csv',index_col=0)
    df_img['tissue'] = [item.split('-Scene')[0] for item in df_img.scene]
    #inconsistient naming in single vs. cyclic
    df_img['rounds'] = [item.replace('(R6)','R5').replace('(R2)','R1').replace('(R3)','R1').replace('(R4)','R1').replace('(R5)','R1') for item in df_img.rounds]
    #add exposure times
    for s_tissue in sorted(set(df_img.tissue)):
        df_exp = pd.read_csv(f'{s_path}/{s_tissue}_ExposureTimes.csv',index_col=0)
        df_tissue = df_img.loc[df_img.tissue==s_tissue,:]
        if s_tissue=='Her2B-K154':
            df_exp = df_exp.append(pd.Series(data={'0':20,'1':50,'2':200,'3':500,'4':500}, name='R11Q_L488.L555.L647.L750_Her2B-K154_2019_2019_01_15__6762-Stitching-01-Scene-06.czi'))
        df_tissue = mpimage.add_exposure(df_tissue, df_exp, type='czi')
        df_img_exp = df_img_exp.append(df_tissue)
    df_img_exp['round_int'] = [float(item.replace('Q','.5').replace('r','.5').split('R')[1]) for item in df_img_exp.rounds]
    df_img_exp.sort_values(['round_int','scene','color'],inplace=True)
    df_img_exp.to_csv(f'{codedir}/Metadata/{s_sample}/{s_sample}_TifImage_ExposureTimes.csv',header=True,index=True)
'''
# 4 combine metadata parameters used for segmentation into a csv 
# this is only necessary for 44290 single vs cyclic and BM-Her2N75 3 TMA replicates
# then copy, save with thresh_* prefix, and enter thresholds for image thresholding analysis  

'''
d_condition = {'44290-112':'cyclic_R1-5_tumor',
 '44290-113':'single_R2_tumor',
 '44290-114':'single_R3_tumor',
 '44290-115':'single_R4_tumor',
 '44290-116':'single_R5_tumor',
 '44294-116':'cyclic_R1-5_normal',
 '44294-117':'single_R2_normal',
 '44294-118':'single_R3_normal',
 '44294-119':'single_R4_normal',
 '44294-120':'single_R5_normal',
 }
df_t = pd.DataFrame()
os.chdir(f'{codedir}/Thresholds')
d_rename = {'HER2':'Her2', 'pHH3':'pH3', 'LamAC':'LaminAC'}
for s_sample, s_condition in d_condition.items():
    print(f'metadata_{s_sample}_RoundsCyclesTable.txt')
    df_tt = pd.read_csv(
        f'metadata_{s_sample}_RoundsCyclesTable.txt',
        delim_whitespace=True,
        header=None,
        names=['marker', 'rounds', 'color', 'minimum', 'max', 'exposure', 'refexp', 'location'],
    )
    df_tt.marker = df_tt.marker.replace(to_replace=d_rename)
    df_tt = df_tt.set_index(f'{s_sample}_' + df_tt.index.astype(str))
    df_tt['slide'] = s_sample
    df_tt['condition'] = d_condition[s_sample]
    df_t = df_t.append(df_tt)
df_t.to_csv('metadata_single_vs_cyclic.csv')

#TMA replicates
ls_sample =  ['BM-Her2N75-15','BM-Her2N75-17','BM-Her2N75-18']
d_rename = {'Her2':'HER2','Lam':'LamAC','pH3':'pHH3'}
df_t = pd.DataFrame()
for s_sample in ls_sample:
    print(f'metadata_{s_sample}_RoundsCyclesTable.txt')
    df_tt = pd.read_csv(
        f'metadata_{s_sample}_RoundsCyclesTable.txt',
        delim_whitespace=True,
        header=None,
        names=['marker', 'rounds', 'color', 'minimum', 'max', 'exposure', 'refexp', 'location'],
    )
    df_tt.marker = df_tt.marker.replace(to_replace=d_rename)
    df_tt = df_tt.set_index(f'{s_sample}_' + df_tt.index.astype(str))
    df_tt['slide'] = s_sample

    df_t = df_t.append(df_tt)
df_t.to_csv('metadata_Jenny_Reps.csv')
'''

# 5 Signal-to-background based on thresholding image

#single versus cyclic
'''
#analyze regions in crop (20191209)
d_crop = {'44290-112':(5000,8800,1300,1800),
 '44290-113':(5886,8111,1300,1800),
 '44290-114':(5010,9242,1300,1800),
 '44290-115':(5975,10490,1300,1800),
 '44290-116':(6336,9174,1300,1800),
 '44294-116':(9547,5459,1300,1800),
 '44294-117':(12180,6092, 1300,1800),
 '44294-118':(10901, 6582, 1300,1800),
 '44294-119':(9853,5557, 1300,1800),
 '44294-120':(10250,5510, 1300,1800),
 }

d_rename = {#'Her2':'HER2', 'pH3':'pHH3', 'LaminAC':'LamAC',
 '(R2)':'R2', '(R3)':'R3', '(R4)':'R4', '(R5)':'R5', '(R6)':'R6'}
df_thresh = pd.read_csv(f'{codedir}/Thresholds/thresh_single_vs_cyclic.csv',index_col=0)

df_thresh['slide_marker'] = df_thresh.slide + '_' + df_thresh.marker

d_process = {
    '44290':f'{codedir}/Images/tiff/44290',
 }
for idx,(s_sample, s_path) in enumerate(d_process.items()):
    df_result = pd.DataFrame()
    os.chdir(s_path)
    df_img = mpimage.parse_org()
    df_img = df_img.replace(d_rename)
    df_img['slide_marker'] = df_img.scene + '_' + df_img.marker
    df_img['img_index'] = df_img.index
    df_thresh = df_img.merge(df_thresh, how='inner', on='slide_marker',suffixes=('','_y'))
    df_thresh.index = df_thresh.img_index
    #results
    df_result = pd.DataFrame()
    dd_result = {}
    for s_marker in sorted(set(df_thresh.marker)):
        
        df_marker = df_thresh[(df_thresh.marker==s_marker) & (df_thresh.rounds !='R6')]
        print(f'{s_marker}  {len(df_marker)}')
        #plot the images
        fig = mpimage.array_img(df_marker.sort_values('condition'),s_xlabel='marker',ls_ylabel=['scene','color'],s_title='condition',tu_array=(2,len(df_marker)//2),tu_fig=(8,8),cmap='inferno')
        fig.savefig(f'{codedir}/Figures/{s_sample}/SinglevsCyclic_TissueLoss_Background_{s_marker}.png')
        df_marker_thresh,d_mask = analyze.thresh_meanint(df_marker.sort_values('condition'),d_crop)
        df_result =df_result.append(df_marker_thresh)
        #plot the mask
        fig, ax = plt.subplots(2,len(df_marker)//2,figsize=(6,8))
        ax=ax.ravel()
        for idx,(s_index, a_mask) in enumerate(d_mask.items()):
            ax[idx].imshow(a_mask)
            ax[idx].set_title(df_thresh.loc[s_index,'condition'])
        fig.savefig(f'{codedir}/Figures/{s_sample}/SinglevsCyclic_Thresholding_{s_marker}.png')

df_result.to_csv(f'{codedir}/Metadata/{s_sample}/SNR_single_vs_cyclic.csv')
'''

#Reproducibility #analyze full TMA core
'''
d_crop ={
    'BM-Her2N75-15-Scene-017':(1160,500,4000,4000),
    'BM-Her2N75-17-Scene-017':(1107,1095,4000,4000),
    'BM-Her2N75-18-Scene-017':(2189,2082,4000,4000),

    'BM-Her2N75-15-Scene-049':(2066,879,4000,4000),
    'BM-Her2N75-17-Scene-049':(2009,1397,4000,4000),
    'BM-Her2N75-18-Scene-049':(1148,699,4000,4000),

    'BM-Her2N75-15-Scene-059':(460,950,4000,4000),
    'BM-Her2N75-17-Scene-059':(449,1353,4000,4000),
    'BM-Her2N75-18-Scene-059':(1400,700,4000,4000),

  }

#load thresholds, rename some columns
df_thresh = pd.read_csv(f'{codedir}/Thresholds/thresh_Jenny_replicates.csv',index_col=0)
df_thresh['scene_num'] = 'Scene-' + df_thresh.scene.astype('str')
df_thresh['scene'] = df_thresh.tissue
df_thresh['minimum'] = df_thresh.threshold
d_process = {
    'BM-Her2N75':f'{codedir}/Images/tiff',
 }
#results
df_result = pd.DataFrame()

for idx,(s_sample, s_path) in enumerate(d_process.items()):
    df_thresh['path'] = f'{s_path}/' + df_thresh.slide + '/' + df_thresh.index
    df_thresh.index = df_thresh.path
    
    for s_marker in sorted(set(df_thresh.marker)):
        
        df_marker = df_thresh[(df_thresh.marker==s_marker)]
        print(f'{s_marker}  {len(df_marker)}')
        #plot the images
        fig = mpimage.array_img(df_marker.sort_values(['scene_num','slide']),s_xlabel='marker',ls_ylabel=['slide','scene_num'],
            s_title='scene',tu_array=(3,len(df_marker)//3),tu_fig=(9,8),cmap='inferno',d_crop=d_crop)
        fig.savefig(f'{codedir}/Figures/{s_sample}/Replicate_TissueLoss_Background_{s_marker}.png')
        #df_marker_thresh,d_mask = analyze.thresh_meanint(df_marker,d_crop)
        #df_result =df_result.append(df_marker_thresh)
        #plot the mask
        #fig, ax = plt.subplots(3,len(df_marker)//3,figsize=(10,10))
        #ax=ax.ravel()
        #for idx,(s_index, a_mask) in enumerate(d_mask.items()):
        #    ax[idx].imshow(a_mask)
        #    ax[idx].set_title(df_thresh.loc[s_index,'scene'])
        #fig.savefig(f'{codedir}/Figures/{s_sample}/Replicate_Thresholding_{s_marker}.png')
        
#df_result.to_csv(f'{codedir}/Metadata/{s_sample}/SNR_jenny_replicates.csv')
'''

#6 fluorescence intensity (quenching experiments
'''
d_process = {
    '4165NPanc':f'{codedir}/Metadata/4165NPanc',
    'NPPan4165':f'{codedir}/Metadata/NPPan4165',
    'B1640':f'{codedir}/Metadata/B1640',
    }

df_result=pd.DataFrame()
for idx,(s_sample, s_path) in enumerate(d_process.items()):
    df_img = pd.read_csv(f'{codedir}/Metadata/{s_sample}/{s_sample}_TifImage_ExposureTimes.csv',index_col=0)
    os.chdir(f'{codedir}/Images/tiff/{s_sample}')
    if s_sample=='NPPan4165':
        df_dapi_r1 = df_img[(df_img.rounds=='R0') & (df_img.color=='c1')]
    else:
        df_dapi_r1 = df_img[(df_img.rounds=='R1') & (df_img.color=='c1')]
    for s_index in df_dapi_r1.index:
        a_dapi = skimage.io.imread(s_index)
        a_dapi_thresh = a_dapi>500
        s_tissue = df_dapi_r1.loc[s_index,'tissue']
        df_tissue = df_img[df_img.tissue==s_tissue]
        df_tissue_result = analyze.mask_meanint(df_tissue, a_mask=a_dapi_thresh)
        #save dapi mask
        skimage.io.imsave(f'{codedir}/Figures/{s_sample}/area_of_measurement_{s_tissue}.png',arr=(a_dapi_thresh.astype('uint8')*255))
        df_result=df_result.append(df_tissue_result)
    df_result.to_csv(f'{codedir}/Metadata/{s_sample}/{s_sample}_MeanIntensityMeasurement.csv',header=True, index=True)
'''
os.chdir(codedir)

