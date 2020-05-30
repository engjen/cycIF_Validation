# visualize segmenatation borders as tifs
# date: 2020-05-30
# author: engje, bue
# language Python 3.6
# license: GPL>=v3

#load libraries
import pandas as pd
import os
import json
from mplex_image import analyze, preprocess
import skimage

############# Paths ############
codedir = '/home/groups/graylab_share/OMERO.rdsStore/engje/Data/cycIF_ValidationStudies'
segdir = f'{codedir}/Images/seg/44290-146'
imagedir = f'{codedir}/Images/tiff'
borderdir = f'{imagedir}/44290-146_R0/Borders'
#parameters
s_pos_name = '20200515_BC44290-146_R0'
#s_pos_name = '20200515_BC44290-146_R8Q'
#borderdir = f'{imagedir}/44290-146_R8/Borders'


preprocess.cmif_mkdir([borderdir])

#load positive data

df_data = pd.read_csv(f'{codedir}/Data/{s_pos_name}_ManualPositive.csv',index_col=0)
print(f'Loaded {s_pos_name}_ManualPositive.csv')
df_data['scene'] = [item.split('_')[1] for item in df_data.index]
df_data['slide'] = [item.split('_')[0] for item in df_data.index]

#Calculate borders/cells touching
'''
ls_color = df_data.columns[df_data.dtypes=='bool']
#ls_color = ['PCNA_Nuclei','CK7_Ring','Ki67_Nuclei','HER2_Ring']
for s_sample in ['BC44290-146']:
    df_pos = df_data[df_data.slide==s_sample]
    analyze.make_border(s_sample,df_pos,ls_color,segdir,savedir=borderdir,b_images=True)
'''

os.chdir(codedir)

