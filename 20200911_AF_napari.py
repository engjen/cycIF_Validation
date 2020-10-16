# 2020-05-05 Need to threshold with Napari!
# author engje
# ipython --gui=qt
#%run 20200504_JPTMAs_napari.py
import napari
import os
import skimage
from skimage import io
import numpy as np
import copy
import pandas as pd
import tifffile

#paths
codedir = os.getcwd()
s_slide = '4165NPanc-74-Scene-001'
regdir = f'{codedir}\\Images\\tiff\\4165NPanc\\Cropped'

#os.chdir('..')
from mplex_image import visualize as viz
from mplex_image import process, analyze

#load positive and intensity data
os.chdir(f'{codedir}\\Data')
df_pos = pd.read_csv(f'4165NPanc_AF_norm_leiden.csv',index_col=0)
df_pos = analyze.celltype_to_bool(df_pos, 'leiden')
df_pos.columns = [str(item) for item in df_pos.columns]

#load images
os.chdir(regdir)
s_crop = 'x2000y10000'  #
#s_crop = 'x2000y7000'
#s_crop = 'x6000y3000'
viewer = napari.Viewer()
label_image = viz.load_crops(viewer,s_crop,s_slide)

#show positive results
ls_cell = df_pos.columns.tolist()
# positive in scene
df_pos['tissue'] = [item.split('_')[0] for item in df_pos.index]
df_scene = df_pos[df_pos.tissue==s_slide.split('-Scene')[0]]

for s_cell in ls_cell: 
    label_image_cell = viz.pos_label(viewer,df_scene,label_image,s_cell)

os.chdir(codedir)


