# 2020-05-05 Need to threshold with Napari!
# author engje
# ipython --gui=qt
# %run 20200902_Annotate_napari.pyrun 20200504_JPTMAs_napari.py
import napari
import os
import skimage
from skimage import io
import numpy as np
import copy
import pandas as pd
import tifffile
import json
from mplex_image import visualize as viz
from mplex_image import process

#paths
codedir = 'C:\\Users\\engje\\Documents\\Data\\PipelineExample'
s_slide = 'BC44290-146-Scene-1'
regdir = f'{codedir}\\44290-146_Cropped'


#load positive and intensity data
os.chdir(codedir)
df_mi = pd.read_csv(f'features_BC44290-146_FilteredMeanIntensity_DAPI8Q_DAPI2.csv',index_col=0)
df_pos = pd.read_csv(f'20200911_BC44290-146_PositiveNegative.csv',index_col=0)
df_af =df_pos.loc[:,df_pos.columns.str.contains('FalseNegative') | df_pos.columns.str.contains('AF')]
df_af.rename({'AF':'AF_'},axis=1,inplace=True)

#load images
os.chdir(regdir)
ls_crop = ['x1800y9000'] 
ls_crop = ['x2000y5000']
ls_crop = ['x2000y7000']

viewer = napari.Viewer()
label_image = viz.load_crops(viewer,ls_crop)


#show positive results
ls_cell = df_af.columns.tolist()

for s_cell in ls_cell: 
    label_image_cell = viz.pos_label(viewer,df_af,label_image,s_cell)


#true negatives (R8Q false positives)

d_result = { #['x2000y7000']
    'CD4R':[[1381.18178064, 1167.52534939],
       [ 571.87259166, 14.32131056]],
    'CD68':[[1619.57331007, 1051.80155136],
       [1447.56806507, 1819.63925294],
       [1283.27001693, 1564.30065155],
       [ 369.84348738, 1597.61376619]],
    'CD8R':[[ 951.82290682,   90.3360093 ],
       [1196.48355621, 1855.65557084]],
    'CD31':[[  82.14010145, 1845.35591055]],
    'CD4':[[1563.24442559,   10.30673509]],
    }

d_result = { #['x2000y5000']
    'CD68':[[1407.54191671, 1983.84245846]],
    #'PDPN':[[]],
    #'aSMA':[[]],
    #'PCNA':[[]],
    'ER':[[619.82496889,  40.60634881]],
    #'CD8':[[]],
    'PD1':[[1033.48474961,   31.36584275],
       [1107.61847063,   27.48110627],
       [1325.4254353 , 1828.14697367],
       [1594.8071995 , 1885.21599927]],
    'CD4':[[1249.57500705,   46.35869932],
       [1965.13887629,   14.71219476]],
    }

#drop lamB1 and Vim
d_result = { #['x1800y9000'] 2000X2000 pixels
    #'AR':[[]],
    #'CD4R': [[]],
    'CD8R':[[1097.51591281, 1464.731752],
       [1726.59810008, 1198.07145173],
       [1633.68815278, 1331.44740395],
       [1783.57423501, 1476.80129363]],
    'CD31':[[1073.83394211, 1584.78914659]],
    #'PDPN':[[]],
    #'CK17':[[]],
    #'Ki67':[[]],
    'CD68':[[ 601.54176882, 1898.32020524],
       [ 628.38177528, 1898.95925302],
       [1588.01111769, 1228.85081276],
       [1477.62764205, 1489.64693652]],
    #'aSMA':[[]],
    #'ER':[[]],
    #'HER2':[[]],
    #'PCNA':[[]],
    'CD4':[[924.8928818 , 307.04451231]],
    #'PD1':[[]],
    #'CK7':[[]],
    #'CD44':[[]],
    #'CK14':[[]],'pHH3':[[]],
    'CD8':[[1865.69698867, 1258.66083733],
       [1950.51373498, 1197.15250985]],
    }

with open(f'TrueNegatives_44290-146_scaledAF_{ls_crop[0]}.json','w') as f: 
    json.dump(d_result, f)

d_cell_result = {}
for s_marker, points in d_result.items():
    ls_points = []
    for point in points:
        point = np.array(point).astype('int')
        i_cell = label_image[point[0]][point[1]]
        ls_points.append(str(i_cell))
    d_cell_result.update({s_marker:ls_points})

with open(f'TrueNegatives_CellIDs_44290-146_scaledAF_{ls_crop[0]}.json','w') as f: 
    json.dump(d_cell_result, f)

os.chdir(codedir)