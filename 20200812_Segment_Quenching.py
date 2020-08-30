# image processing for with mlpex_image
# date: 2020-07-21
# author: engje
# language Python 3.7
# license: GPL>=v3

#libraries

import os
import sys
import numpy as np
import pandas as pd
import shutil
import matplotlib.pyplot as plt

# cd /home/groups/graylab_share/OMERO.rdsStore/engje/Data/cycIF_ValidationStudies/cycIF_Validation/Images/tiff/4165NPanc

####  Paths  ####

codedir = os.getcwd()
#czidir = codedir.replace('_Workflow','_Images')
#tiffdir = f'{codedir}/RawImages'
#qcdir = f'{codedir}/QC'
regdir = f'{codedir}/Images/tiff/4165NPanc'
subdir = f'{regdir}/SubtractedRegisteredImages'
segdir = f'{regdir}/Segmentation'

# Start Preprocessing

from mplex_image import preprocess, mpimage, cmif
preprocess.cmif_mkdir([segdir,subdir]) #tiffdir,qcdir,regdir,

ls_sample = [
 '4165NPanc-73',
 '4165NPanc-75',
 '4165NPanc-77',
 #'4165NPanc-74',
 #'4165NPanc-76',
 #'4165NPanc-78',
 #'B1640-T8-3',  'B1640-T8-5', #resegment on eppec tmux 0 (3-cell); tmux 1 (5-nuc); tmux2 5-cell
 ]

#### 2 Tiffs: check/change names ####
'''
for s_sample in ls_sample: 
    os.chdir(f'{subdir}/{s_sample}')
    #Example: change file name and change back
    d_rename = {'R1_blank.blank.blank.blank':'R1_R1c2.R1c3.R1c4.R1c5',
        'R2_blank.blank.blank.blank':'R2_R2c2.R2c3.R2c4.R2c5', 
        'R3_blank.blank.blank.blank':'R3_R3c2.R3c3.R3c4.R3c5',
        'R4_blank.blank.blank.blank':'R4_R4c2.R4c3.R4c4.R4c5',
        'R5_blank.blank.blank.blank':'R5_R5c2.R5c3.R5c4.R5c5',
        'R6_blank.blank.blank.blank':'R6_R6c2.R6c3.R6c4.R6c5'}
    preprocess.dchange_fname(d_rename)#,b_test=False)
    #sort and count images 
    df_img = mpimage.parse_org() 
    #cmif.count_images(df_img)
'''

### 7 move images and prep for seg ###
'''
for s_sample in ls_sample:
     preprocess.cmif_mkdir([f'{subdir}/{s_sample}'])

'''
#### 8 cell pose segmentation & feat extraction
'''
from mplex_image import segment
import time

nuc_diam = 30
cell_diam = 30 

s_seg_markers= "['Ecad']" #"['CK7-Vimentin']"#"['blank']" ###
s_type = 'nuclei'#'cell'#


print(f'Predicting {s_type}')
for s_sample in ls_sample:
    preprocess.cmif_mkdir([f'{segdir}/{s_sample}Cellpose_Segmentation'])
    os.chdir(f'{subdir}')
    for s_file in os.listdir():
        if s_file.find(s_sample) > -1:
            os.chdir(f'{subdir}/{s_file}')
            print(f'Processing {s_file}')
            df_img = segment.parse_org()
            for s_scene in sorted(set(df_img.scene)):
                s_slide_scene= f'{s_sample}-Scene-{s_scene}'
                s_find = df_img[(df_img.rounds=='R1') & (df_img.color=='c1') & (df_img.scene==s_scene)].index[0]
                segment.cellpose_segment_job(s_sample,s_slide_scene,s_find,f'{codedir}/Images/tiff/4165NPanc',nuc_diam,cell_diam,s_type,s_seg_markers,b_match=True) 
                os.chdir(f'{segdir}/{s_sample}Cellpose_Segmentation')
                os.system(f'sbatch cellpose_{s_type}_{s_slide_scene}.sh')
                time.sleep(12)
                print('Next')
'''
#### Extract Cellpose Features ####### 

from mplex_image import features

nuc_diam = 30
cell_diam = 30 
ls_seg_markers = ['Ecad']#['blank']#

for s_sample in ls_sample: 
    df_sample, df_thresh = features.extract_cellpose_features(s_sample, segdir, subdir, ls_seg_markers, nuc_diam, cell_diam)
    df_sample.to_csv(f'{segdir}/features_{s_sample}_MeanIntensity_Centroid_Shape.csv')
    df_thresh.to_csv(f'{segdir}thresh_{s_sample}_ThresholdLi.csv')

#### 9 filter cellpose data 
'''
nuc_diam = 20
cell_diam = 25 

ls_seg_markers = ['Ecad']

d_job ={#'JE-TMA-41':['Scene 03_','Scene 04_','Scene 05_','Scene 09_'],
    'SMT101Bx2-3':['Scene 001','Scene 002'],
    'BC44290-146':['Scene 01_'],
  }

#parameters
ls_marker_cyto = [
        'CK14_cytoplasm','CK5_cytoplasm','CK17_cytoplasm',
        'CK19_cytoplasm','CK7_cytoplasm',#'pAKT_cytoplasm',#'CK8_cytoplasm','EGFR_cytoplasm',
        'Ecad_cytoplasm','HER2_cytoplasm']

ls_color = ['DAPI1_nuclei', 'DAPI2_nuclei', 'DAPI3_nuclei', 'DAPI4_nuclei',
 'DAPI5_nuclei', 'DAPI6_nuclei','DAPI7_nuclei', 'DAPI8_nuclei', #'DAPI9_nuclei', 'DAPI10_nuclei','DAPI11_nuclei'
 ]

ls_filter =  ['DAPI3_nuclei','DAPI6_nuclei']

#filter
df_mi_full = process.load_cellpose_df(ls_sample, segdir)
df_xy = process.filter_cellpose_xy(df_mi_full)
#df_mi_full = df_mi_full.loc[:,[len(item.split('_')) < 3 for item in df_mi_full.columns]] #hack
#df_mi_full = df_mi_full.loc[:,((~df_mi_full.columns.str.contains('_cell')) & (~df_mi_full.columns.str.contains('DAPI11_')))]#hack because of cell
#df_img_all = process.load_li_thresh(ls_sample, segdir)
#df_mi_filled = process.fill_cellpose_nas(df_img_all,df_mi_full,ls_marker_cyto,s_thresh='Ecad',ls_celline=[],ls_shrunk = ['CD44','Vim'])
#df_mi,df_mi_nas = process.filter_cellpose_df(df_mi_filled)
#df_pos_auto,d_thresh_record = process.auto_threshold(df_mi,df_img_all)
#process.autothresh_scatterplots(df_pos_auto,d_thresh_record,df_mi,df_xy,ls_color)
#df_mi_filter = process.filter_dapi_cellpose(df_pos_auto,ls_color,df_mi,ls_filter,df_img_all)

# Output the data 
#print('filtered columns:')
#print(set(df_mi_filter.columns))
s_sample  = codedir.split('cmIF_')[1]
#df_mi_filter.to_csv(f'{segdir}/features_{s_sample}_FilteredMeanIntensity_{"_".join([item.split("_")[0] for item in ls_filter])}.csv')
df_xy.to_csv(f'{segdir}/features_{s_sample}_CentroidXY.csv')

#combine labels for cells touching analysis

for s_sample, ls_scene in d_job.items():
    se_neg = df_mi_full[df_mi_full.slide == s_sample].Ecad_negative
    labels,combine,dd_result = features.combine_labels(s_sample,ls_scene, segdir, subdir, ls_seg_markers, nuc_diam, cell_diam, se_neg)
'''

#edge mask
'''
i_pixel = 154

for s_sample, ls_scene in d_job.items(): #{s_sample}Cellpose_Segmentation
    df_img = pd.read_csv(f'{segdir}/thresh_{s_sample}_ThresholdLi.csv',index_col=0)
    for s_scene in ls_scene:
        s_index = df_img[(df_img.index.str.contains(s_scene.replace(' ','-'))) & (df_img.index.str.contains('R1_')) & (df_img.index.str.contains('_c1_'))].index[0]
        img_dapi = io.imread(f'{subdir}/{s_sample}/{s_index}')
        mask = img_dapi > 350 #df_img.loc[s_index,'threshold_li']
        mask_filled = morphology.remove_small_holes(mask, 20000)
        border_mask, __, __,distances = features.mask_border(mask_filled,type='inner',pixel_distance = i_pixel)
        img = np.zeros(border_mask.shape,dtype='uint8')
        img[border_mask] = 255
        io.imsave(f"{segdir}/TissueEdgeMask{i_pixel}_{s_sample}_scene{s_scene.split(' ')[1]}.png", img)
        

#edge ROI (Cellpose)

os.chdir(f'{segdir}')

for s_sample, ls_scene in d_job.items():
    s_sample_data  = codedir.split('cmIF_')[1]
    #load xy
    df_sample = pd.DataFrame()
    df_xy = pd.read_csv(f'features_{s_sample_data}_CentroidXY.csv',index_col=0)
    df_xy['cells'] = [int(item.split('cell')[1]) for item in df_xy.index]
    #load masks
    for s_scene in ls_scene:
        mask = io.imread(f"{segdir}/TissueEdgeMask{i_pixel}_{s_sample}_scene{s_scene.split(' ')[1]}.png")
        mask_gray = mask == 255
        labels = io.imread(f'{segdir}/{s_sample}Cellpose_Segmentation/{s_scene} nuclei{nuc_diam} - Nuclei Segmentation Basins.tif')
        edge = features.mask_labels(mask_gray,labels)
        df_scene = df_xy[df_xy.index.str.contains(f'{s_sample}_scene{s_scene.split(" ")[1].split("_")[0]}')]
        #old (use if you have coordinates, not labels)
        #mask_gray = mask#[:,:,0]
        #contour = skimage.measure.find_contours(mask_gray,0)
        #coords = skimage.measure.approximate_polygon(contour[0], tolerance=5)
        #fig,ax=plt.subplots()
        #ax.imshow(mask_gray)
        #ax.plot(coords[:, 1], coords[:, 0], '-r', linewidth=2)
        #fig.savefig(f'TissueEdgeMask_{s_sample}_Scene-{s_scene}_polygon.png')
        #x = np.array(df_scene.DAPI_X.astype('int').values)
        #y = np.array(df_scene.DAPI_Y.astype('int').values)
        #points = np.array((y,x)).T
        #mask = skimage.measure.points_in_poly(points, coords)
        #works
        es_cells = set(edge.astype('int')).intersection(set(df_scene.cells))
        df_edge = df_scene[df_scene.cells.isin(es_cells)]
        fig,ax=plt.subplots()
        ax.imshow(mask_gray)
        ax.scatter(df_edge.DAPI_X,df_edge.DAPI_Y,s=1)
        fig.savefig(f'TissueEdgeMask_{s_sample}_Scene-{s_scene}_cells.png')
        df_sample = df_sample.append(df_edge)

    df_sample.to_csv(f'features_{s_sample_data}_EdgeCells{i_pixel}pixels_CentroidXY.csv')
'''


#10 PNGs
'''
#20200423 make multicolor png for figure 3 (800x800 pixel)
d_overlay = {'R1':['CD20','CD8','CD4','CK19'],
     'R2':[ 'PCNA','HER2','ER','CD45'],
     'R3':['pHH3', 'CK14', 'CD44', 'CK5'],
     'R4':[ 'Vim', 'CK7', 'CD31', 'LamAC',],
     'R5':['aSMA', 'CD68', 'Ki67', 'Ecad'],
     'R6':['CK17','PDPN','CD31','CD3'],
     'R7':['CD20P','CD8R','CD4R','CK5P'],
     'R8':['LamB1','AR','ColIV','ColI'],
     'subtype':[ 'PCNA','HER2','ER','Ki67'],
     'diff':['Ecad', 'CK14', 'CD44', 'CK5'],
     'immune':['aSMA','CD8','CD4','CK19'],
     'stromal':['aSMA','Vim','CD44','CK19'],
     }

d_crop ={
    #'JE-TMA-41-Scene-03':(2100,2300),
    #'JE-TMA-41-Scene-04':(1916,2603),
    #'JE-TMA-41-Scene-05':(3800,2100),
    #'JE-TMA-41-Scene-09':(2000,2000),
    #'JE-TMA-41-Scene-08':(2000,2000),
    'JE-TMA-41-Scene-07':(1500,1500),
    #'JE-TMA-41-Scene-06':(2000,2000),
    #'JE-TMA-41-Scene-10':(2200,2000),
    #'JE-TMA-41-Scene-11':(2200,2000),
    #'BC44290-146-Scene-01':(2287,9202),#(3140,3946),
    #'SMT101Bx2-3-Scene-001':(1380,2756),
    #'SMT101Bx2-3-Scene-002':(4002,3986),
  }

#registered (not subtracted)

es_bright = {'pHH3','CD68','CK14','CK5','CK17'} #tissue
#es_bright = {'pHH3','CD68','CK14','CK5','CK17','CD3','PDPN','CD45'} #SMT tissue

#es_bright = {'CD20','CD8','CD4','pHH3','CD45','aSMA',
# 'PDPN','CD31','CD3','CD68','CK14','CK5','CK17','ColIV','ColI'} #cell line
'''

'''
for s_scene,item in d_crop.items():
    preprocess.cmif_mkdir([f'{qcdir}/{s_scene}'])
    s_path = f'{regdir}/{s_scene}'
    os.chdir(s_path)
    df_img = mpimage.parse_org()
    df_img['path'] = [f'{s_path}/{item}' for item in df_img.index]
    for s_scene in sorted(set(df_img.scene)):
        df_dapi = df_img[(df_img.color=='c1')&(df_img.scene==s_scene)]
        df_scene = df_img[(df_img.color!='c1') & (df_img.scene==s_scene)]
        for s_round in ['stromal','diff','subtype','R1','R2','R3','R4','R5','R6','R7','R8']:
            #df_round = df_scene[df_scene.rounds==s_round]
            df_round = df_scene[df_scene.marker.isin(d_overlay[s_round])]
            #df_dapi_round = df_dapi[(df_dapi.rounds==s_round)]
            df_dapi_round = df_dapi[(df_dapi.rounds=='R2')]
            high_thresh=0.999
            d_overlay_round = {s_round:d_overlay[s_round]}
            d_result = mpimage.multicolor_png(df_round,df_dapi_round,s_scene=s_scene,d_overlay=d_overlay_round,d_crop=d_crop,es_dim={'nada'},es_bright=es_bright,low_thresh=2000,high_thresh=high_thresh)
            for key, tu_result in d_result.items():
                skimage.io.imsave(f'{qcdir}/{s_scene}/ColorArray_{s_scene}_{key}_{".".join(tu_result[0])}.png',tu_result[1])
'''
#AF subtracted
'''
for s_sample in ls_sample:
    preprocess.cmif_mkdir([f'{qcdir}/{s_sample}'])
    s_path = f'{subdir}/{s_sample}'
    os.chdir(s_path)
    df_img = mpimage.parse_org()
    df_img['path'] = [f'{s_path}/{item}' for item in df_img.index]
    for s_scene in sorted(set(df_img.scene)):
        df_dapi = df_img[(df_img.color=='c1')&(df_img.scene==s_scene)]
        df_scene = df_img[(df_img.color!='c1') & (df_img.scene==s_scene)]
        for s_round in ['R1','R2','R3','R4','R5','R6','R7','R8']:
            df_round = df_scene[df_scene.rounds==s_round]
            df_dapi_round = df_dapi[(df_dapi.rounds==s_round)]
            high_thresh=0.999
            d_overlay_round = {s_round:d_overlay[s_round]}
            d_result = mpimage.multicolor_png(df_round,df_dapi_round,s_scene=s_scene,d_overlay=d_overlay_round,d_crop=d_crop,es_dim={'nada'},es_bright=es_bright,low_thresh=2000,high_thresh=high_thresh)
            for key, tu_result in d_result.items():
                skimage.io.imsave(f'{qcdir}/{s_sample}/ColorArray_{s_scene}_{key}_{".".join(tu_result[0])}.png',tu_result[1])
'''
os.chdir(codedir)

