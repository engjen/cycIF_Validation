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
qcdir = f'{regdir}/QC'
cropdir = f'{regdir}/Cropped'

# Start Preprocessing

from mplex_image import preprocess, mpimage, cmif
preprocess.cmif_mkdir([segdir,subdir,qcdir,cropdir]) #tiffdir,,regdir,

ls_sample = [
 #'4165NPanc-73',
 #'4165NPanc-75',
 #'4165NPanc-77',
 #'4165NPanc-74',
 #'4165NPanc-76',
 '4165NPanc-78',
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

s_seg_markers= "['CK7-Vimentin']"#"['Ecad']" #"['blank']" ###
s_type = 'cell'#'nuclei'#


print(f'Predicting {s_type}')
for s_sample in ls_sample:
    preprocess.cmif_mkdir([f'{segdir}/{s_sample}Cellpose_Segmentation'])
    #os.chdir(f'{subdir}')
    #os.chdir(f'{subdir}')
    #for s_file in os.listdir():
    #    if s_file.find(s_sample) > -1:
    #        os.chdir(f'{subdir}/{s_file}')
    #        print(f'Processing {s_file}')
    os.chdir(regdir)
    df_img = segment.parse_org()
    for s_scene in sorted(set(df_img.scene)):
                s_slide_scene= f'{s_sample}-Scene-{s_scene}'
                s_find = df_img[(df_img.rounds=='R1') & (df_img.color=='c1') & (df_img.scene==s_scene)].index[0]
                segment.cellpose_segment_job(s_sample,s_slide_scene,
                 s_find,f'{segdir}/{s_sample}Cellpose_Segmentation',
                 f'{regdir}',nuc_diam,cell_diam,
                 s_type,s_seg_markers,b_long=True)#,b_match=True) #b_long=True)                os.chdir(f'{segdir}/{s_sample}Cellpose_Segmentation')
                os.system(f'sbatch cellpose_{s_type}_{s_slide_scene}.sh')
                time.sleep(9)
                print('Next')
'''

#### Extract Cellpose Features ####### 
'''
from mplex_image import features

nuc_diam = 30
cell_diam = 30 
ls_seg_markers = ['Ecad']#['blank']#

for s_sample in ls_sample: 
    df_sample, df_thresh = features.extract_cellpose_features(s_sample, segdir, subdir, ls_seg_markers, nuc_diam, cell_diam)
    df_sample.to_csv(f'{segdir}/features_{s_sample}_MeanIntensity_Centroid_Shape.csv')
    df_thresh.to_csv(f'{segdir}thresh_{s_sample}_ThresholdLi.csv')
'''
#### 9 filter cellpose data 
'''
from mplex_image import process
#parameters
ls_filter = ['DAPI2_nuclei','DAPI6_nuclei']
ls_marker_cyto = []
ls_custom = []
s_thresh= 'Ecad'

for s_sample in ls_sample:
    #filtering
    os.chdir(segdir)
    df_img_all = process.load_li([s_sample])
    #dapi filter too high
    #df_img_all.loc[df_img_all[df_img_all.marker.str.contains('DAPI')].index,'threshold_li'] = 500
    df_mi_full = process.load_cellpose_df([s_sample], segdir)
    df_xy = process.filter_cellpose_xy(df_mi_full)
    df_mi_filled = process.fill_cellpose_nas(df_img_all,df_mi_full,ls_marker_cyto,s_thresh=s_thresh,ls_celline=[],ls_shrunk = [],qcdir=qcdir)
    df_mi = process.filter_loc_cellpose(df_mi_filled, ls_marker_cyto, ls_custom)
    df_pos_auto,d_thresh_record = process.auto_threshold(df_mi,df_img_all)
    ls_color = df_mi.columns[df_mi.columns.str.contains('DAPI')].tolist()
    process.positive_scatterplots(df_pos_auto,d_thresh_record,df_xy,ls_color,qcdir) #+ [f'{s_thresh}_cytoplasm']
    df_mi_filter = process.filter_dapi_cellpose(df_pos_auto,ls_color,df_mi,ls_filter,df_img_all,qcdir)
    df_mi_filter.to_csv(f'{segdir}/features_{s_sample}_FilteredMeanIntensity_{"_".join([item.split("_")[0] for item in ls_filter])}.csv')
    df_xy.to_csv(f'{segdir}/features_{s_sample}_CentroidXY.csv')

#Expand nuclei without matching cell labels for cells' touching analysis
#for s_sample in ls_sample:
#    se_neg = df_mi_full[df_mi_full.slide==s_sample].Ecad_negative
#    labels,combine,dd_result = features.combine_labels(s_sample, segdir, subdir, ls_seg_markers, nuc_diam, cell_diam, se_neg)

'''

d_crop = {
 '4165NPanc-74-Scene-001': (6000,3000),
 '4165NPanc-76-Scene-001': (6000,3000),
 '4165NPanc-78-Scene-001': (6000,3000),
 }
 
d_crop = {
 #'4165NPanc-74-Scene-001': (2000,10000),
 #'4165NPanc-76-Scene-001': (2000,10000),
 '4165NPanc-78-Scene-001': (2000,10000),
 }
tu_dim=(2000,2000)

#10 -2 ome-tiff
import tifffile
#ome-tiff parameters
s_dapi = 'DAPI1'
d_combos = {#'AF':{'R0c2','R0c3','R0c4','R0c5'}, 
        'Background':{'R1c2','R2c2','R3c2','R4c2','R5c2','R1c3','R1c4'},
    }
for s_sample in ls_sample:
    os.chdir(f'{subdir}/{s_sample}')
    df_img = mpimage.parse_org()
    for s_index in df_img.index:
        s_marker =  df_img.loc[s_index,'marker']
        if s_marker == 'DAPI':
            s_marker = s_marker + f'{df_img.loc[s_index,"rounds"].split("R")[1]}'
    dd_result = mpimage.overlay_crop(d_combos,d_crop,df_img,s_dapi,tu_dim)
    for s_crop, d_result in dd_result.items():
        for s_type, (ls_marker, array) in d_result.items():
            new_array = array[np.newaxis,np.newaxis,:]
            s_xml =  mpimage.gen_xml(new_array, ls_marker)
            with tifffile.TiffWriter(f'{cropdir}/{s_crop}_{s_type}.ome.tif') as tif:
                tif.save(new_array,  photometric = "minisblack", description=s_xml, metadata = None)

#10-3 crop basins to match cropped overlays

#cmif.crop_labels(d_crop,tu_dim,segdir,cropdir,s_find='exp5_CellSegmentationBasins')
#cmif.crop_labels(d_crop,tu_dim,segdir,cropdir,s_find='Nuclei Segmentation Basins')

os.chdir(codedir)

