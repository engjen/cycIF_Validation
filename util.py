import os
import skimage
from skimage import io, segmentation, morphology, measure
import numpy as np
import tifffile
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import pandas as pd

import patsy
import sys
import numpy.linalg as la
import sys

import re
import shutil
import itertools
from itertools import chain
import json
os.chdir('/home/groups/graylab_share/OMERO.rdsStore/engje/Data/mplex_image')
import ometiff

#functions

def parse_org(s_end = "ORG.tif",s_start='R',type='reg'):
    """
    This function will parse images following koei's naming convention
    Example: Registered-R1_PCNA.CD8.PD1.CK19_Her2B-K157-Scene-002_c1_ORG.tif
    The output is a dataframe with image filename in index
    And rounds, color, imagetype, scene (/tissue), and marker in the columns
    type= 'reg' or 'raw'
    """

    ls_file = []
    for file in os.listdir():
    #find all filenames ending in s_end
        if file.endswith(s_end):
            if file.find(s_start)==0:
                ls_file = ls_file + [file]
    lls_name = [item.split('_') for item in ls_file]
    df_img = pd.DataFrame(index=ls_file)
    if type == 'raw':
        lls_scene = [item.split('-Scene-') for item in ls_file]
    elif type== 'noscenes':
        ls_scene = ['Scene-001'] * len(ls_file)
    if type == 'raw':
        df_img['rounds'] = [item[0] for item in lls_name]
    elif type== 'noscenes':
        df_img['rounds'] = [item[0] for item in lls_name]
    else:
        df_img['rounds'] = [item[0].split('Registered-')[1] for item in lls_name]
    df_img['color'] = [item[-2] for item in lls_name]
    df_img['imagetype'] = [item[-1].split('.tif')[0] for item in lls_name]
    if type == 'raw':
        df_img['slide'] = [item[2] for item in lls_name]
        try:
            df_img['scene'] = [item[1].split('_')[0] for item in lls_scene]
        except IndexError:
            print(f"{set([item[0] for item in lls_scene])}")
    elif type == 'noscenes':
        df_img['slide'] = [item[2] for item in lls_name]
        df_img['scene'] = ls_scene
    else:
        df_img['scene'] = [item[2] for item in lls_name]
    df_img['round_ord'] = [re.sub('Q','.5', item) for item in df_img.rounds] 
    df_img['round_ord'] = [float(re.sub('[^0-9.]','', item)) for item in df_img.round_ord]
    df_img = df_img.sort_values(['round_ord','rounds','color'])
    for idx, s_round in enumerate(df_img.rounds.unique()):
        df_img.loc[df_img.rounds==s_round, 'round_num'] = idx
    #parse file name for biomarker
    for s_index in df_img.index:
        #print(s_index)
        s_color = df_img.loc[s_index,'color']
        if s_color == 'c1':
            s_marker = 'DAPI'
        elif s_color == 'c2':
            s_marker = s_index.split('_')[1].split('.')[0]
        elif s_color == 'c3':
            s_marker = s_index.split('_')[1].split('.')[1]
        elif s_color == 'c4':
            s_marker = s_index.split('_')[1].split('.')[2]
        elif s_color == 'c5':
            s_marker = s_index.split('_')[1].split('.')[3]
        #these are only included in sardana shading corrected images
        elif s_color == 'c6':
            s_marker = s_index.split('_')[1].split('.')[2]
        elif s_color == 'c7':
            s_marker = s_index.split('_')[1].split('.')[3]
        else: print('Error')
        df_img.loc[s_index,'marker'] = s_marker

    return(df_img) #,lls_name)

def array_roi_if(df_img,df_dapi,s_label='rounds',s_title='Title',tu_crop=(0,0,100,100),tu_array=(2,4),tu_fig=(10,20),tu_rescale=(0,0),i_expnorm=0,i_micron_per_pixel=.325):
    """
    create a grid of images
    df_img = dataframe of images with columns having image attributes
        and index with image names
    df_dapi = like df_img, but with the matching dapi images
    s_label= attribute to label axes
    s_title = x axis title
    tu_crop = (upper left corner x,  y , xlength, yheight)
    tu_array = subplot array dimensions
    tu_fig = size of figue
    tu_rescale= range of rescaling
    i_expnorm = normalize to an exposure time (requires 'exposure' column in dataframe
    """
    cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', [(0,0,0),(0,1,0)], N=256, gamma=1.0)
    fig, ax = plt.subplots(tu_array[0],tu_array[1],figsize=tu_fig,sharey=True, squeeze=False) #
    ax = ax.ravel()
    for ax_num, s_index in enumerate(df_img.index):
        s_col_label = df_img.loc[s_index,s_label]
        #load image, copr, rescale
        a_image=io.imread(s_index)
        a_dapi = io.imread((df_dapi).index[0])# & (df_dapi.rounds=='R1')
        a_crop = a_image[(tu_crop[1]):(tu_crop[1]+tu_crop[3]),(tu_crop[0]):(tu_crop[0]+tu_crop[2])]
        a_crop_dapi = a_dapi[(tu_crop[1]):(tu_crop[1]+tu_crop[3]),(tu_crop[0]):(tu_crop[0]+tu_crop[2])]
        #a_crop_dapi = (a_crop_dapi/255).astype('int')
        if i_expnorm > 0:
            a_crop = a_crop/df_img.loc[s_index,'exposure']*i_expnorm
        if tu_rescale==(0,0):
            a_rescale = skimage.exposure.rescale_intensity(a_crop,in_range=(np.quantile(a_crop,0.03),1.5*np.quantile(a_crop,0.998)),out_range=(0, 255))
            tu_max = (np.quantile(a_crop,0.03),1.5*np.quantile(a_crop,0.998))
        else:
            #print(f'original {a_crop.min()},{a_crop.max()}')
            #print(f'rescale to {tu_rescale}')
            a_rescale = skimage.exposure.rescale_intensity(a_crop,in_range = tu_rescale,out_range=(0,255))
            tu_max=tu_rescale
        a_rescale_dapi = skimage.exposure.rescale_intensity(a_crop_dapi,in_range = (np.quantile(a_crop_dapi,0.03),2*np.quantile(a_crop_dapi,0.99)),out_range=(0,255)) 
        a_rescale_dapi = a_rescale_dapi.astype(np.uint8)
        a_rescale = a_rescale.astype(np.uint8)
        #2 color png
        zdh = np.dstack((np.zeros_like(a_rescale), a_rescale, a_rescale_dapi))
        ax[ax_num].imshow(zdh)
        ax[ax_num].set_title('')
        ax[ax_num].set_ylabel('')
        ax[ax_num].set_xlabel(s_col_label,fontsize = 'x-large')
        if tu_rescale == (0,0):
            if len(ax)>1:
                ax[ax_num].set_xlabel(f'{s_col_label} ({int(np.quantile(a_crop,0.03))} - {int(1.5*np.quantile(a_crop,0.998))})')
        ax[ax_num].set_xticklabels('')
    #pixel to micron (apply after ax is returned)
    #ax[0].set_yticklabels([str(int(re.sub(u"\u2212", "-", item.get_text()))*i_micron_per_pixel) for item in ax[0].get_yticklabels(minor=False)])
    plt.suptitle(s_title,y=0.95,size = 'xx-large',weight='bold')
    plt.subplots_adjust(wspace=.05, hspace=.05)
    # Now adding the colorbar
    norm = mpl.colors.Normalize(vmin=tu_max[0],vmax=tu_max[1])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    if len(ax) == 1:
        cbaxes = fig.add_axes([.88, 0.125, 0.02, 0.75]) #[left, bottom, width, height]
        plt.colorbar(sm, cax=cbaxes)#,format=ticker.FuncFormatter(fmt))
        plt.figtext(0.47,0.03,s_label.replace('_',' '),fontsize = 'x-large', weight='bold')
    elif tu_rescale != (0,0):
        cbaxes = fig.add_axes([.91, 0.15, 0.015, 0.7]) #[left, bottom, width, height]
        plt.colorbar(sm, cax=cbaxes)#,format=ticker.FuncFormatter(fmt))
        plt.figtext(0.42,0.03,s_label.replace('_',' '),fontsize = 'x-large', weight='bold')
    else:
        print("Different ranges - can't use colorbar") 
        plt.figtext(0.43,0.03,s_label.replace('_',' '),fontsize = 'x-large', weight='bold')

    return(fig,ax) 

def array_img(df_img,s_xlabel='color',ls_ylabel=['rounds','exposure'],s_title='marker',tu_array=(2,4),tu_fig=(10,20),cmap='gray',d_crop={}):
    """
    create a grid of images
    df_img = dataframe of images with columns having image attributes
        and index with image names
    s_xlabel = coumns of grid
    ls_ylabel = y label 
    s_title= title

    """
     
    fig, ax = plt.subplots(tu_array[0],tu_array[1],figsize=tu_fig)
    ax = ax.ravel()
    for ax_num, s_index in enumerate(df_img.index):
        s_row_label = f'{df_img.loc[s_index,ls_ylabel[0]]}\n {df_img.loc[s_index,ls_ylabel[1]]}'
        s_col_label = df_img.loc[s_index,s_xlabel]
        a_image=skimage.io.imread(s_index)
        s_label_img = df_img.loc[s_index,s_title]
        a_rescale = skimage.exposure.rescale_intensity(a_image,in_range=(0,1.5*np.quantile(a_image,0.98)))
        if len(d_crop)!= 0:
            tu_crop = d_crop[df_img.loc[s_index,'scene']]
            a_rescale = a_rescale[(tu_crop[1]):(tu_crop[1]+tu_crop[3]),(tu_crop[0]):(tu_crop[0]+tu_crop[2])]
        ax[ax_num].imshow(a_rescale,cmap=cmap)
        ax[ax_num].set_title(s_label_img)
        ax[ax_num].set_ylabel(s_row_label)
        ax[ax_num].set_xlabel(f'{s_col_label}\n 0 - {int(1.5*np.quantile(a_image,0.98))}')
    plt.tight_layout()
    return(fig)

def thresh_meanint(df_thresh,d_crop={},s_thresh='minimum',):
    """
    threshold, and output positive and negative mean intensity and array
    df_thresh = dataframe of images with columns having image attributes
        and index with image names, column with threshold values
    d_crop = image scene and crop coordinates

    """
    d_mask = {}
    for idx, s_index in enumerate(df_thresh.index):
        #load image, crop, thresh
        a_image = skimage.io.imread(s_index)
        if len(d_crop) != 0:
            tu_crop = d_crop[df_thresh.loc[s_index,'scene']]
            a_image = a_image[(tu_crop[1]):(tu_crop[1]+tu_crop[3]),(tu_crop[0]):(tu_crop[0]+tu_crop[2])]
        i_min = df_thresh.loc[s_index,s_thresh]
        a_mask = a_image > i_min
        print(f'mean positive intensity = {np.mean(a_image[a_mask])}')
        df_thresh.loc[s_index,'meanpos'] = np.mean(a_image[a_mask])
        b_mask = a_image < i_min
        print(f'mean negative intensity = {np.mean(a_image[b_mask])}')
        df_thresh.loc[s_index,'meanneg'] = np.mean(a_image[b_mask])
        d_mask.update({s_index:a_mask})
    return(df_thresh,d_mask)

def quartiles(regionmask, intensity):
    return np.percentile(intensity[regionmask], q=(5,25, 50, 75,95))

def thresh_erode(df_thresh,d_crop={},s_thresh='minimum',k=10):
    """
    threshold, erode around pixels above threshold to obtain background,
    and output foreground and background intensity 
    df_thresh = dataframe of images with columns having image attributes
        and index with image names, column with threshold values
    d_crop = image scene and crop coordinates

    """
    d_mask = {}
    df_all = pd.DataFrame()
    for idx, s_index in enumerate(df_thresh.index):
        #load image, crop, thresh
        a_image = skimage.io.imread(s_index)
        if len(d_crop) != 0:
            tu_crop = d_crop[df_thresh.loc[s_index,'scene']]
            a_image = a_image[(tu_crop[1]):(tu_crop[1]+tu_crop[3]),(tu_crop[0]):(tu_crop[0]+tu_crop[2])]
        # generate foreground  
        fg = a_image > df_thresh.loc[s_index,s_thresh]
        fg = skimage.morphology.remove_small_objects(fg,min_size=120) #remove flecks
        #generate background
        bg = fg==0
        bg = morphology.binary_erosion(bg, morphology.disk(30)) #30 pixels - 10 um (one cell diam)
        bg = morphology.remove_small_objects(bg,min_size=1000)  #a couple cells in size
        #superpixels to label
        suppix_b = segmentation.slic(np.ones(bg.shape), n_segments=k,start_label=1)   
        suppix_b[~bg] = 0
        suppix_f = segmentation.slic(np.ones(bg.shape), n_segments=k,start_label=1) 
        suppix_f[~fg] = 0
        #measure
        props_f = measure.regionprops_table(label_image=suppix_f,intensity_image=a_image,properties=('label','mean_intensity','centroid'),extra_properties=(quartiles,))
        props_b = measure.regionprops_table(label_image=suppix_b,intensity_image=a_image,properties=('label','mean_intensity','centroid'),extra_properties=(quartiles,))
        df = pd.DataFrame(props_b).merge(pd.DataFrame(props_f),on='label',suffixes=('_bg','_fg'))
        df['marker'] = df_thresh.loc[s_index,'marker']
        df['filename'] = s_index
        print(f'SBR {df_thresh.loc[s_index,"marker"]} {(df.mean_intensity_fg/df.mean_intensity_bg).mean()}')
        df_all = df_all.append(df)        
        d_mask.update({s_index:(suppix_f,suppix_b)})
    return(df_all,d_mask)

def overlay_crop(d_combos,d_crop,df_img,s_dapi,tu_dim=(1000,1000),b_8bit=True): 
    """
    output custon multi page tiffs according to dictionary, with s_dapi as channel 1 in each overlay
    BUG with 53BP1
    d_crop : {slide_scene : (x,y) coord
    tu_dim = (width, height)
    d_combos = {'Immune':{'CD45', 'PD1', 'CD8', 'CD4', 'CD68', 'FoxP3','GRNZB','CD20','CD3'},
    'Stromal':{'Vim', 'aSMA', 'PDPN', 'CD31', 'ColIV','ColI'},
    'Differentiation':{'CK19', 'CK7','CK5', 'CK14', 'CK17','CK8'},
    'Tumor':{'HER2', 'Ecad', 'ER', 'PgR','Ki67','PCNA'},
    'Proliferation':{'EGFR','CD44','AR','pHH3','pRB'}, 
    'Functional':{'pS6RP','H3K27','H3K4','cPARP','gH2AX','pAKT','pERK'},
    'Lamins':{'LamB1','LamAC', 'LamB2'}}
    """
    dd_result = {}
    for s_index in df_img.index:
        s_marker =  df_img.loc[s_index,'marker']
        if s_marker == 'DAPI':
            s_marker = s_marker + f'{df_img.loc[s_index,"rounds"].split("R")[1]}'
        df_img.loc[s_index,'marker'] = s_marker
    #now make overlays
    for s_scene, xy_cropcoor in d_crop.items():
        d_result = {}
        print(f'Processing {s_scene}')
        df_slide = df_img[df_img.slide_scene==s_scene]
        s_image_round = df_slide[df_slide.marker==s_dapi].index[0]
        if len(df_slide[df_slide.marker==s_dapi.split('_')[0]].index) == 0:
            print('Error: dapi not found')
        elif len(df_slide[df_slide.marker==s_dapi.split('_')[0]].index) > 1:
            print('Error: too many dapi images found')
        else:
            print(s_image_round)
        #exclude any missing biomarkers
        es_all = set(df_slide.marker)
        #iterate over overlay combinations
        for s_type, es_combos in d_combos.items():
            d_overlay = {}
            es_combos_shared = es_combos.intersection(es_all)
            for idx, s_combo in enumerate(sorted(es_combos_shared)):
                s_filename = (df_slide[df_slide.marker==s_combo]).index[0]
                if len((df_slide[df_slide.marker==s_combo]).index) == 0:
                    print(f'Error: {s_combo} not found')
                elif len((df_slide[df_slide.marker==s_combo]).index) > 1:
                    print(f'\n Warning {s_combo}: too many marker images found, used {s_filename}')
                else:
                    print(f'{s_combo}: {s_filename}')
                d_overlay.update({s_combo:s_filename})
            #d_overlay.update({s_dapi:s_image_round})
            a_dapi = io.imread(s_image_round)
            #crop 
            a_crop = a_dapi[(xy_cropcoor[1]):(xy_cropcoor[1]+tu_dim[1]),(xy_cropcoor[0]):(xy_cropcoor[0]+tu_dim[0])]
            a_overlay = np.zeros((len(d_overlay) + 1,a_crop.shape[0],a_crop.shape[1]),dtype=np.uint8)
            if a_crop.dtype == 'uint16':
                if b_8bit:
                    a_crop = (a_crop/256).astype(np.uint8)
                else:
                    a_rescale = skimage.exposure.rescale_intensity(a_crop,in_range=(0,1.5*np.quantile(a_crop,0.9999)))
                    a_crop = (a_rescale/256).astype(np.uint8)
                    print(f'rescale intensity')
            a_overlay[0,:,:] = a_crop
            ls_biomarker_all = [s_dapi]
            for i, s_color in enumerate(sorted(d_overlay.keys())):
                s_overlay= d_overlay[s_color]
                ls_biomarker_all.append(s_color)
                a_channel = io.imread(s_overlay)
                #crop 
                a_crop = a_channel[(xy_cropcoor[1]):(xy_cropcoor[1]+tu_dim[1]),(xy_cropcoor[0]):(xy_cropcoor[0]+tu_dim[0])]
                if a_crop.dtype == 'uint16':
                    if b_8bit:
                        a_crop = (a_crop/256).astype(np.uint8)
                    else:
                        a_rescale = skimage.exposure.rescale_intensity(a_crop,in_range=(0,1.5*np.quantile(a_crop,0.9999)))
                        a_crop = (a_rescale/256).astype(np.uint8)
                        print(f'rescale intensity')
                a_overlay[i + 1,:,:] = a_crop
            d_result.update({s_type:(ls_biomarker_all,a_overlay)})
        dd_result.update({f'{s_scene}_x{xy_cropcoor[0]}y{xy_cropcoor[1]}':d_result})
        return(dd_result)
    
def cropped_ometiff(dd_result,cropdir):
    for s_crop, d_result in dd_result.items():
        for s_type, (ls_marker, array) in d_result.items():
            print(f'Generating multi-page ome-tiff {[item for item in ls_marker]}')
            new_array = array[np.newaxis,np.newaxis,:]
            s_xml =  ometiff.gen_xml(new_array, ls_marker)
            with tifffile.TiffWriter(f'{cropdir}/{s_crop}_{s_type}.ome.tif') as tif:
                tif.save(new_array,  photometric = "minisblack", description=s_xml, metadata = None)
                
def multicolor_png(df_img,df_dapi,s_scene,d_overlay,d_crop,es_dim={'CD8','FoxP3','ER','AR'},es_bright={'Ki67','pHH3'},low_thresh=4000,high_thresh=0.999):
    '''
    create RGB image with Dapi plus four - 6 channels
    '''

    d_result = {}
    #print(s_scene)
    tu_crop = d_crop[s_scene]
    df_slide = df_img[df_img.scene == s_scene]
    x=tu_crop[1]
    y=tu_crop[0]
    img_dapi = skimage.io.imread(df_dapi[df_dapi.scene==s_scene].path[0])
    a_crop = img_dapi[x:x+800,y:y+800]
    a_rescale_dapi = skimage.exposure.rescale_intensity(a_crop,in_range=(np.quantile(img_dapi,0.2),1.5*np.quantile(img_dapi,high_thresh)),out_range=(0, 255))
    if 1.5*np.quantile(img_dapi,high_thresh) < low_thresh:
                a_rescale_dapi = skimage.exposure.rescale_intensity(a_crop,in_range=(low_thresh/2,low_thresh),out_range=(0, 255))
    elif len(es_dim.intersection(set(['DAPI'])))==1:
                new_thresh = float(str(high_thresh)[:-2])
                a_rescale_dapi = skimage.exposure.rescale_intensity(a_crop,in_range=(np.quantile(img_dapi,0.2),1.5*np.quantile(img_dapi,new_thresh)),out_range=(0, 255))
    elif len(es_bright.intersection(set(['DAPI'])))==1:
                a_rescale_dapi = skimage.exposure.rescale_intensity(a_crop,in_range=(np.quantile(img_dapi,0.2),1.5*np.quantile(img_dapi,float(str(high_thresh) + '99'))),out_range=(0, 255))

    #RGB
    for s_type, ls_marker in d_overlay.items():
        #print(s_type)
        zdh = np.dstack((np.zeros_like(a_rescale_dapi), np.zeros_like(a_rescale_dapi),a_rescale_dapi))
        for idx, s_marker in enumerate(ls_marker):
            #print(s_marker)
            s_index = df_slide[df_slide.marker == s_marker].index[0]
            img = skimage.io.imread(df_slide.loc[s_index,'path'])
            a_crop = img[x:x+800,y:y+800]
            in_range = (np.quantile(a_crop,0.2),1.5*np.quantile(a_crop,high_thresh))
            a_rescale = skimage.exposure.rescale_intensity(a_crop,in_range=in_range,out_range=(0, 255))
            if 1.5*np.quantile(a_crop,high_thresh) < low_thresh:
                #print('low thresh')
                in_range=(low_thresh/2,low_thresh)
                a_rescale = skimage.exposure.rescale_intensity(a_crop,in_range=in_range,out_range=(0, 255))
            elif len(es_dim.intersection(set([s_marker])))==1:
                #print('dim')
                new_thresh = float(str(high_thresh)[:-2])
                in_range=(np.quantile(a_crop,0.2),1.5*np.quantile(a_crop,new_thresh))
                a_rescale = skimage.exposure.rescale_intensity(a_crop,in_range=in_range,out_range=(0, 255))
            elif len(es_bright.intersection(set([s_marker])))==1:
                #print('bright')
                in_range=(np.quantile(a_crop,0.2),1.5*np.quantile(a_crop,float(str(high_thresh) + '99')))
                a_rescale = skimage.exposure.rescale_intensity(a_crop,in_range=in_range,out_range=(0, 255))

            #print(f'low {int(in_range[0])} high {int(in_range[1])}')
            if idx == 0:
                zdh = zdh + np.dstack((np.zeros_like(a_rescale), a_rescale,np.zeros_like(a_rescale)))

            elif idx == 1:
                zdh = zdh + np.dstack((a_rescale, a_rescale,np.zeros_like(a_rescale)))

            elif idx == 2:
                zdh = zdh + np.dstack((a_rescale, np.zeros_like(a_rescale),np.zeros_like(a_rescale) ))

            elif idx == 3:
                zdh = zdh + np.dstack((np.zeros_like(a_rescale), a_rescale, a_rescale))
        #print(zdh.min())
        zdh = zdh.clip(0,255)
        zdh = zdh.astype('uint8')
        #print(zdh.max())
        d_result.update({s_type:(ls_marker,zdh)})
    return(d_result)


def roi_if_border(df_img,df_dapi,df_border,s_label='rounds',s_title='Title',tu_crop=(0,0,100,100),tu_array=(2,4),tu_fig=(10,20),tu_rescale=(0,0),i_expnorm=0,i_micron_per_pixel=.325):
    """
    create a grid of images
    df_img = dataframe of images with columns having image attributes
        and index with image names
    df_dapi = like df_img, but with the matching dapi images
    df_border: index is border image file name
    s_label= attribute to label axes
    s_title = x axis title
    tu_crop = (upper left corner x,  y , xlength, yheight)
    tu_array = subplot array dimensions
    tu_fig = size of figue
    tu_rescale= 
    i_expnorm = 
    """
    cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', [(0,0,0),(0,1,0)], N=256, gamma=1.0)
    fig, ax = plt.subplots(tu_array[0],tu_array[1],figsize=tu_fig,sharey=True, squeeze=False) #
    ax = ax.ravel()
    for ax_num, s_index in enumerate(df_img.index):
        s_col_label = df_img.loc[s_index,s_label]
        #load image, copr, rescale
        a_image=skimage.io.imread(s_index)
        a_dapi = skimage.io.imread((df_dapi).index[0])# & (df_dapi.rounds=='R1')
        a_crop = a_image[(tu_crop[1]):(tu_crop[1]+tu_crop[3]),(tu_crop[0]):(tu_crop[0]+tu_crop[2])]
        a_crop_dapi = a_dapi[(tu_crop[1]):(tu_crop[1]+tu_crop[3]),(tu_crop[0]):(tu_crop[0]+tu_crop[2])]
        #a_crop_dapi = (a_crop_dapi/255).astype('int')
        if i_expnorm > 0:
            a_crop = a_crop/df_img.loc[s_index,'exposure']*i_expnorm
        if tu_rescale==(0,0):
            a_rescale = skimage.exposure.rescale_intensity(a_crop,in_range=(np.quantile(a_crop,0.03),1.5*np.quantile(a_crop,0.998)),out_range=(0, 255))
            tu_max = (np.quantile(a_crop,0.03),1.5*np.quantile(a_crop,0.998))
        else:
            print(f'original {a_crop.min()},{a_crop.max()}')
            print(f'rescale to {tu_rescale}')
            a_rescale = skimage.exposure.rescale_intensity(a_crop,in_range = tu_rescale,out_range=(0,255))
            tu_max=tu_rescale
        a_rescale_dapi = skimage.exposure.rescale_intensity(a_crop_dapi,in_range = (np.quantile(a_crop_dapi,0.03),2*np.quantile(a_crop_dapi,0.99)),out_range=(0,255)) 
        a_rescale_dapi = a_rescale_dapi.astype(np.uint8)
        a_rescale = a_rescale.astype(np.uint8)
        #white border
        s_border_index = df_border[df_border.marker==(df_img.loc[s_index,'marker'])].index[0]
        a_border = skimage.io.imread(s_border_index)
        a_crop_border = a_border[(tu_crop[1]):(tu_crop[1]+tu_crop[3]),(tu_crop[0]):(tu_crop[0]+tu_crop[2])]
        mask = a_crop_border > 250
        #2 color png
        zdh = np.dstack((np.zeros_like(a_rescale), a_rescale, a_rescale_dapi))
        zdh[mask] = 255
        #zdh = zdh.clip(0,255)
        #zdh = zdh.astype('uint8')
        ax[ax_num].imshow(zdh)
        ax[ax_num].set_title('')
        ax[ax_num].set_ylabel('')
        ax[ax_num].set_xlabel(s_col_label,fontsize = 'x-large')
        if tu_rescale == (0,0):
            if len(ax)>1:
                ax[ax_num].set_xlabel(f'{s_col_label} ({int(np.quantile(a_crop,0.03))} - {int(1.5*np.quantile(a_crop,0.998))})')
        ax[ax_num].set_xticklabels('')
    #pixel to micron (apply after ax is returned)
    #ax[0].set_yticklabels([str(int(re.sub(u"\u2212", "-", item.get_text()))*i_micron_per_pixel) for item in ax[0].get_yticklabels(minor=False)])
    plt.suptitle(s_title,y=0.93,size = 'xx-large',weight='bold')
    plt.subplots_adjust(wspace=.05, hspace=.05)
    # Now adding the colorbar
    norm = mpl.colors.Normalize(vmin=tu_max[0],vmax=tu_max[1])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    if len(ax) == 1:
        cbaxes = fig.add_axes([.88, 0.125, 0.02, 0.75]) #[left, bottom, width, height]
        plt.colorbar(sm, cax = cbaxes)
        plt.figtext(0.47,0.03,s_label.replace('_',' '),fontsize = 'x-large', weight='bold')
    elif tu_rescale != (0,0):
        cbaxes = fig.add_axes([.92, 0.175, 0.02, 0.64]) #[left, bottom, width, height]
        plt.colorbar(sm, cax = cbaxes)
        plt.figtext(0.42,0.03,s_label.replace('_',' '),fontsize = 'x-large', weight='bold')
    else:
        print("Different ranges - can't use colorbar") 
        plt.figtext(0.43,0.03,s_label.replace('_',' '),fontsize = 'x-large', weight='bold')

    return(fig,ax,a_crop_border) 

def add_exposure(df_img,df_t,type='roundcycles'):
    """
    df_img = dataframe of images with columns [ 'color', 'exposure', 'marker','sub_image','sub_exposure']
            and index with image names
    df_t = metadata with dataframe with ['marker','exposure']
    """
    if type == 'roundscycles':
        for s_index in df_img.index:
            s_marker = df_img.loc[s_index,'marker']
            #look up exposure time for marker in metadata
            df_t_image = df_t[(df_t.marker==s_marker)]
            if len(df_t_image) > 0:
                i_exposure = df_t_image.iloc[0].loc['exposure']
                df_img.loc[s_index,'exposure'] = i_exposure
            else:
                print(f'{s_marker} has no recorded exposure time')
    elif type == 'czi':
    #add exposure
        df_t['rounds'] = [item.split('_')[0] for item in df_t.index]
        #df_t['tissue'] = [item.split('_')[2].split('-Scene')[0] for item in df_t.index] #not cool with stiched 
        for s_index in df_img.index:
            s_tissue = df_img.loc[s_index,'scene'].split('-Scene')[0]
            s_color = str(int(df_img.loc[s_index,'color'].split('c')[1])-1)
            s_round = df_img.loc[s_index,'rounds']
            print(s_index)
            df_img.loc[s_index,'exposure'] = df_t[(df_t.index.str.contains(s_tissue)) & (df_t.rounds==s_round)].loc[:,s_color][0]

    return(df_img)

def cmif_mkdir(ls_dir):
    '''
    check if directories existe. if not, make them
    '''
    for s_dir in ls_dir:
        if not os.path.exists(s_dir):
            os.makedirs(s_dir)

#normalization functions
def aprior(gamma_hat):
    m = gamma_hat.mean()
    s2 = gamma_hat.var()
    return (2 * s2 +m**2) / s2

def bprior(gamma_hat):
    m = gamma_hat.mean()
    s2 = gamma_hat.var()
    return (m*s2+m**3)/s2

def it_sol(sdat, g_hat, d_hat, g_bar, t2, a, b, conv=0.0001):
    n = (1 - np.isnan(sdat)).sum(axis=1)
    g_old = g_hat.copy()
    d_old = d_hat.copy()

    change = 1
    count = 0
    while change > conv:
        #print g_hat.shape, g_bar.shape, t2.shape
        g_new = postmean(g_hat, g_bar, n, d_old, t2)
        sum2 = ((sdat - np.dot(g_new.values.reshape((g_new.shape[0], 1)), np.ones((1, sdat.shape[1])))) ** 2).sum(axis=1)
        d_new = postvar(sum2, n, a, b)
       
        change = max((abs(g_new - g_old) / g_old).max(), (abs(d_new - d_old) / d_old).max())
        g_old = g_new #.copy()
        d_old = d_new #.copy()
        count = count + 1
    adjust = (g_new, d_new)
    return adjust 

def postmean(g_hat, g_bar, n, d_star, t2):
    return (t2*n*g_hat+d_star * g_bar) / (t2*n+d_star)

def postvar(sum2, n, a, b):
    return (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)

def design_mat(mod, numerical_covariates, batch_levels):
    # require levels to make sure they are in the same order as we use in the
    # rest of the script.
    design = patsy.dmatrix("~ 0 + C(batch, levels=%s)" % str(batch_levels),
                                                  mod, return_type="dataframe")

    mod = mod.drop(["batch"], axis=1)
    numerical_covariates = list(numerical_covariates)
    sys.stderr.write("found %i batches\n" % design.shape[1])
    other_cols = [c for i, c in enumerate(mod.columns)
                  if not i in numerical_covariates]
    factor_matrix = mod[other_cols]
    design = pd.concat((design, factor_matrix), axis=1)
    if numerical_covariates is not None:
        sys.stderr.write("found %i numerical covariates...\n"
                            % len(numerical_covariates))
        for i, nC in enumerate(numerical_covariates):
            cname = mod.columns[nC]
            sys.stderr.write("\t{0}\n".format(cname))
            design[cname] = mod[mod.columns[nC]]
    sys.stderr.write("found %i categorical variables:" % len(other_cols))
    sys.stderr.write("\t" + ", ".join(other_cols) + '\n')
    return design

def combat(data, batch, model=None, numerical_covariates=None):
    """Correct for batch effects in a dataset
    Parameters
    ----------
    data : pandas.DataFrame
        A (n_features, n_samples) dataframe of the expression or methylation
        data to batch correct
    batch : pandas.Series
        A column corresponding to the batches in the data, with index same as
        the columns that appear in ``data``
    model : patsy.design_info.DesignMatrix, optional
        A model matrix describing metadata on the samples which could be
        causing batch effects. If not provided, then will attempt to coarsely
        correct just from the information provided in ``batch``
    numerical_covariates : list-like
        List of covariates in the model which are numerical, rather than
        categorical
    Returns
    -------
    corrected : pandas.DataFrame
        A (n_features, n_samples) dataframe of the batch-corrected data
    """
    if isinstance(numerical_covariates, str):
        numerical_covariates = [numerical_covariates]
    if numerical_covariates is None:
        numerical_covariates = []

    if model is not None and isinstance(model, pd.DataFrame):
        model["batch"] = list(batch)
    else:
        model = pd.DataFrame({'batch': batch})

    batch_items = model.groupby("batch").groups.items()
    batch_levels = [k for k, v in batch_items]
    batch_info = [v for k, v in batch_items]
    n_batch = len(batch_info)
    n_batches = np.array([len(v) for v in batch_info])
    n_array = float(sum(n_batches))

    # drop intercept
    drop_cols = [cname for cname, inter in  ((model == 1).all()).iteritems() if inter == True]
    drop_idxs = [list(model.columns).index(cdrop) for cdrop in drop_cols]
    model = model[[c for c in model.columns if not c in drop_cols]]
    numerical_covariates = [list(model.columns).index(c) if isinstance(c, str) else c
            for c in numerical_covariates if not c in drop_cols]

    design = design_mat(model, numerical_covariates, batch_levels)

    sys.stderr.write("Standardizing Data across genes.\n")
    #error shapes (3,7200) and (26,7200) not aligned: 7200 (dim 1) != 26 (dim 0)
    B_hat = np.dot(np.dot(la.inv(np.dot(design.T, design)), design.T), data.T) #data.T
    grand_mean = np.dot((n_batches / n_array).T, B_hat[:n_batch,:])
    var_pooled = np.dot(((data - np.dot(design, B_hat).T)**2), np.ones((int(n_array), 1)) / int(n_array))

    stand_mean = np.dot(grand_mean.T.reshape((len(grand_mean), 1)), np.ones((1, int(n_array))))
    tmp = np.array(design.copy())
    tmp[:,:n_batch] = 0
    stand_mean  += np.dot(tmp, B_hat).T

    s_data = ((data - stand_mean) / np.dot(np.sqrt(var_pooled), np.ones((1, int(n_array)))))

    sys.stderr.write("Fitting L/S model and finding priors\n")
    batch_design = design[design.columns[:n_batch]]
    gamma_hat = np.dot(np.dot(la.inv(np.dot(batch_design.T, batch_design)), batch_design.T), s_data.T)

    delta_hat = []

    for i, batch_idxs in enumerate(batch_info):
        #batches = [list(model.columns).index(b) for b in batches]
        delta_hat.append(s_data[batch_idxs].var(axis=1))

    gamma_bar = gamma_hat.mean(axis=1) 
    t2 = gamma_hat.var(axis=1)
   

    a_prior = list(map(aprior, delta_hat))
    b_prior = list(map(bprior, delta_hat))

    sys.stderr.write("Finding parametric adjustments\n")
    gamma_star, delta_star = [], []
    for i, batch_idxs in enumerate(batch_info):
        #print '18 20 22 28 29 31 32 33 35 40 46'
        #print batch_info[batch_id]

        temp = it_sol(s_data[batch_idxs], gamma_hat[i],
                     delta_hat[i], gamma_bar[i], t2[i], a_prior[i], b_prior[i])

        gamma_star.append(temp[0])
        delta_star.append(temp[1])

    sys.stdout.write("Adjusting data\n")
    bayesdata = s_data
    gamma_star = np.array(gamma_star)
    delta_star = np.array(delta_star)


    for j, batch_idxs in enumerate(batch_info):

        dsq = np.sqrt(delta_star[j,:])
        dsq = dsq.reshape((len(dsq), 1))
        denom =  np.dot(dsq, np.ones((1, n_batches[j])))
        numer = np.array(bayesdata[batch_idxs] - np.dot(batch_design.loc[batch_idxs], gamma_star).T)

        bayesdata[batch_idxs] = numer / denom
   
    vpsq = np.sqrt(var_pooled).reshape((len(var_pooled), 1))
    bayesdata = bayesdata * np.dot(vpsq, np.ones((1, int(n_array)))) + stand_mean
 
    return bayesdata

#adapted from https://github.com/brentp/combat.py/blob/master/combat.py


def combat_fit(data, batch, model=None, numerical_covariates=None):
    """Correct for batch effects in a dataset
    Parameters
    ----------
    data : pandas.DataFrame
        A (n_features, n_samples) dataframe of the expression or methylation
        data to batch correct
    batch : pandas.Series
        A column corresponding to the batches in the data, with index same as
        the columns that appear in ``data``
    model : patsy.design_info.DesignMatrix, optional
        A model matrix describing metadata on the samples which could be
        causing batch effects. If not provided, then will attempt to coarsely
        correct just from the information provided in ``batch``
    numerical_covariates : list-like
        List of covariates in the model which are numerical, rather than
        categorical
    Returns
    -------
    gamma_star : centering parameters from combat fitting
    delta_star : scaling parameters from combat fitting
    stand_mean: pooled mean of batches
    var_pooled: pooled variance of batches
    """
    if isinstance(numerical_covariates, str):
        numerical_covariates = [numerical_covariates]
    if numerical_covariates is None:
        numerical_covariates = []

    if model is not None and isinstance(model, pd.DataFrame):
        model["batch"] = list(batch)
    else:
        model = pd.DataFrame({'batch': batch})

    batch_items = model.groupby("batch").groups.items()
    batch_levels = [k for k, v in batch_items]
    batch_info = [v for k, v in batch_items]
    n_batch = len(batch_info)
    n_batches = np.array([len(v) for v in batch_info])
    n_array = float(sum(n_batches))

    # drop intercept
    drop_cols = [cname for cname, inter in  ((model == 1).all()).iteritems() if inter == True]
    drop_idxs = [list(model.columns).index(cdrop) for cdrop in drop_cols]
    model = model[[c for c in model.columns if not c in drop_cols]]
    numerical_covariates = [list(model.columns).index(c) if isinstance(c, str) else c
            for c in numerical_covariates if not c in drop_cols]

    design = design_mat(model, numerical_covariates, batch_levels)

    sys.stderr.write("Standardizing Data across genes.\n")
    B_hat = np.dot(np.dot(la.inv(np.dot(design.T, design)), design.T), data.T) 
    grand_mean = np.dot((n_batches / n_array).T, B_hat[:n_batch,:])
    var_pooled = np.dot(((data - np.dot(design, B_hat).T)**2), np.ones((int(n_array), 1)) / int(n_array))

    stand_mean = np.dot(grand_mean.T.reshape((len(grand_mean), 1)), np.ones((1, int(n_array))))
    tmp = np.array(design.copy())
    tmp[:,:n_batch] = 0
    stand_mean  += np.dot(tmp, B_hat).T

    s_data = ((data - stand_mean) / np.dot(np.sqrt(var_pooled), np.ones((1, int(n_array)))))

    sys.stderr.write("Fitting L/S model and finding priors\n")
    batch_design = design[design.columns[:n_batch]]
    gamma_hat = np.dot(np.dot(la.inv(np.dot(batch_design.T, batch_design)), batch_design.T), s_data.T)

    delta_hat = []

    for i, batch_idxs in enumerate(batch_info):
        delta_hat.append(s_data[batch_idxs].var(axis=1))

    gamma_bar = gamma_hat.mean(axis=1) 
    t2 = gamma_hat.var(axis=1)


    a_prior = list(map(aprior, delta_hat))
    b_prior = list(map(bprior, delta_hat))

    sys.stderr.write("Finding parametric adjustments\n")
    gamma_star, delta_star = [], []
    for i, batch_idxs in enumerate(batch_info):
        temp = it_sol(s_data[batch_idxs], gamma_hat[i],
                     delta_hat[i], gamma_bar[i], t2[i], a_prior[i], b_prior[i])

        gamma_star.append(temp[0])
        delta_star.append(temp[1])
    #just retrun one stand_mean array
    stand_mean = stand_mean[:,0]
    return(gamma_star, delta_star, stand_mean, var_pooled)
        
def combat_transform(data, batch, gamma_star, delta_star, stand_mean, var_pooled,model=None, numerical_covariates=None):
    """Correct for batch effects in a dataset
    Parameters
    ----------
    data : pandas.DataFrame
        A (n_features, n_samples) dataframe of the expression or methylation
        data to batch correct
    batch : pandas.Series
        A column corresponding to the batches in the data, with index same as
        the columns that appear in ``data``
    gamma_star : centering parameters from combat fitting
    delta_star : scaling parameters from combat fitting
    stand_mean: pooled mean of batches
    var_pooled: pooled variance of batches
    model : patsy.design_info.DesignMatrix, optional
        A model matrix describing metadata on the samples which could be
        causing batch effects. If not provided, then will attempt to coarsely
        correct just from the information provided in ``batch``
    numerical_covariates : list-like
        List of covariates in the model which are numerical, rather than
        categorical
    Returns
    -------
    corrected : pandas.DataFrame
        A (n_features, n_samples) dataframe of the batch-corrected data
    """
    #get design
    if isinstance(numerical_covariates, str):
        numerical_covariates = [numerical_covariates]
    if numerical_covariates is None:
        numerical_covariates = []

    if model is not None and isinstance(model, pd.DataFrame):
        model["batch"] = list(batch)
    else:
        model = pd.DataFrame({'batch': batch})
    batch_items = model.groupby("batch").groups.items()
    batch_levels = [k for k, v in batch_items]
    batch_info = [v for k, v in batch_items]
    n_batch = len(batch_info)
    n_batches = np.array([len(v) for v in batch_info])
    n_array = float(sum(n_batches))
    # drop intercept
    drop_cols = [cname for cname, inter in  ((model == 1).all()).iteritems() if inter == True]
    drop_idxs = [list(model.columns).index(cdrop) for cdrop in drop_cols]
    model = model[[c for c in model.columns if not c in drop_cols]]
    numerical_covariates = [list(model.columns).index(c) if isinstance(c, str) else c
            for c in numerical_covariates if not c in drop_cols]
    design = design_mat(model, numerical_covariates, batch_levels)
    #standardize
    sys.stderr.write("Standardizing Data across genes.\n")

    #reshape stand mean
    stand_mean = np.dot(stand_mean.T.reshape((len(stand_mean), 1)), np.ones((1, int(data.shape[1]))))
    s_data = ((data - stand_mean) / np.dot(np.sqrt(var_pooled), np.ones((1, int(n_array)))))
    batch_design = design[design.columns[:n_batch]]
    # adjust data
    sys.stdout.write("Adjusting data\n")
    bayesdata = s_data
    gamma_star = np.array(gamma_star)
    delta_star = np.array(delta_star)
    #for each batch
    for j, batch_idxs in enumerate(batch_info):

        dsq = np.sqrt(delta_star[j,:])
        dsq = dsq.reshape((len(dsq), 1))
        denom =  np.dot(dsq, np.ones((1, n_batches[j]))) #divide by sqrt delta_star
        numer = np.array(bayesdata[batch_idxs] - np.dot(batch_design.loc[batch_idxs], gamma_star).T) #subtract gamma_star

        bayesdata[batch_idxs] = numer / denom
    #multiply by square root of variance and add mean
    vpsq = np.sqrt(var_pooled).reshape((len(var_pooled), 1))
    bayesdata = bayesdata * np.dot(vpsq, np.ones((1, int(n_array)))) + stand_mean
    return bayesdata


def combat_fit_old(data, batch, model=None, numerical_covariates=None):
    """Correct for batch effects in a dataset
    Parameters
    ----------
    data : pandas.DataFrame
        A (n_features, n_samples) dataframe of the expression or methylation
        data to batch correct
    batch : pandas.Series
        A column corresponding to the batches in the data, with index same as
        the columns that appear in ``data``
    model : patsy.design_info.DesignMatrix, optional
        A model matrix describing metadata on the samples which could be
        causing batch effects. If not provided, then will attempt to coarsely
        correct just from the information provided in ``batch``
    numerical_covariates : list-like
        List of covariates in the model which are numerical, rather than
        categorical
    Returns
    -------
    gamma_star : centering parameters from combat fitting
    delta_star : scaling parameters from combat fitting
    """
    if isinstance(numerical_covariates, str):
        numerical_covariates = [numerical_covariates]
    if numerical_covariates is None:
        numerical_covariates = []

    if model is not None and isinstance(model, pd.DataFrame):
        model["batch"] = list(batch)
    else:
        model = pd.DataFrame({'batch': batch})

    batch_items = model.groupby("batch").groups.items()
    batch_levels = [k for k, v in batch_items]
    batch_info = [v for k, v in batch_items]
    n_batch = len(batch_info)
    n_batches = np.array([len(v) for v in batch_info])
    n_array = float(sum(n_batches))

    # drop intercept
    drop_cols = [cname for cname, inter in  ((model == 1).all()).iteritems() if inter == True]
    drop_idxs = [list(model.columns).index(cdrop) for cdrop in drop_cols]
    model = model[[c for c in model.columns if not c in drop_cols]]
    numerical_covariates = [list(model.columns).index(c) if isinstance(c, str) else c
            for c in numerical_covariates if not c in drop_cols]

    design = design_mat(model, numerical_covariates, batch_levels)

    sys.stderr.write("Standardizing Data across genes.\n")
    B_hat = np.dot(np.dot(la.inv(np.dot(design.T, design)), design.T), data.T) 
    grand_mean = np.dot((n_batches / n_array).T, B_hat[:n_batch,:])
    var_pooled = np.dot(((data - np.dot(design, B_hat).T)**2), np.ones((int(n_array), 1)) / int(n_array))

    stand_mean = np.dot(grand_mean.T.reshape((len(grand_mean), 1)), np.ones((1, int(n_array))))
    tmp = np.array(design.copy())
    tmp[:,:n_batch] = 0
    stand_mean  += np.dot(tmp, B_hat).T

    s_data = ((data - stand_mean) / np.dot(np.sqrt(var_pooled), np.ones((1, int(n_array)))))

    sys.stderr.write("Fitting L/S model and finding priors\n")
    batch_design = design[design.columns[:n_batch]]
    gamma_hat = np.dot(np.dot(la.inv(np.dot(batch_design.T, batch_design)), batch_design.T), s_data.T)

    delta_hat = []

    for i, batch_idxs in enumerate(batch_info):
        delta_hat.append(s_data[batch_idxs].var(axis=1))

    gamma_bar = gamma_hat.mean(axis=1) 
    t2 = gamma_hat.var(axis=1)
   

    a_prior = list(map(aprior, delta_hat))
    b_prior = list(map(bprior, delta_hat))

    sys.stderr.write("Finding parametric adjustments\n")
    gamma_star, delta_star = [], []
    for i, batch_idxs in enumerate(batch_info):
        temp = it_sol(s_data[batch_idxs], gamma_hat[i],
                     delta_hat[i], gamma_bar[i], t2[i], a_prior[i], b_prior[i])

        gamma_star.append(temp[0])
        delta_star.append(temp[1])
    return(gamma_star, delta_star)
        
def combat_transform_old(data, batch, gamma_star, delta_star,model=None, numerical_covariates=None):
    """Correct for batch effects in a dataset
    Parameters
    ----------
    data : pandas.DataFrame
        A (n_features, n_samples) dataframe of the expression or methylation
        data to batch correct
    batch : pandas.Series
        A column corresponding to the batches in the data, with index same as
        the columns that appear in ``data``
    gamma_star : centering parameters from combat fitting
    delta_star : scaling parameters from combat fitting
    model : patsy.design_info.DesignMatrix, optional
        A model matrix describing metadata on the samples which could be
        causing batch effects. If not provided, then will attempt to coarsely
        correct just from the information provided in ``batch``
    numerical_covariates : list-like
        List of covariates in the model which are numerical, rather than
        categorical
    Returns
    -------
    corrected : pandas.DataFrame
        A (n_features, n_samples) dataframe of the batch-corrected data
    """
    #get design
    if isinstance(numerical_covariates, str):
        numerical_covariates = [numerical_covariates]
    if numerical_covariates is None:
        numerical_covariates = []

    if model is not None and isinstance(model, pd.DataFrame):
        model["batch"] = list(batch)
    else:
        model = pd.DataFrame({'batch': batch})
    batch_items = model.groupby("batch").groups.items()
    batch_levels = [k for k, v in batch_items]
    batch_info = [v for k, v in batch_items]
    n_batch = len(batch_info)
    n_batches = np.array([len(v) for v in batch_info])
    n_array = float(sum(n_batches))
    # drop intercept
    drop_cols = [cname for cname, inter in  ((model == 1).all()).iteritems() if inter == True]
    drop_idxs = [list(model.columns).index(cdrop) for cdrop in drop_cols]
    model = model[[c for c in model.columns if not c in drop_cols]]
    numerical_covariates = [list(model.columns).index(c) if isinstance(c, str) else c
            for c in numerical_covariates if not c in drop_cols]
    design = design_mat(model, numerical_covariates, batch_levels)
    #standardize
    sys.stderr.write("Standardizing Data across genes.\n")
    B_hat = np.dot(np.dot(la.inv(np.dot(design.T, design)), design.T), data.T) 
    grand_mean = np.dot((n_batches / n_array).T, B_hat[:n_batch,:])
    var_pooled = np.dot(((data - np.dot(design, B_hat).T)**2), np.ones((int(n_array), 1)) / int(n_array))

    stand_mean = np.dot(grand_mean.T.reshape((len(grand_mean), 1)), np.ones((1, int(n_array))))
    tmp = np.array(design.copy())
    tmp[:,:n_batch] = 0
    stand_mean  += np.dot(tmp, B_hat).T
    s_data = ((data - stand_mean) / np.dot(np.sqrt(var_pooled), np.ones((1, int(n_array)))))
    batch_design = design[design.columns[:n_batch]]
    # adjust data
    sys.stdout.write("Adjusting data\n")
    bayesdata = s_data
    gamma_star = np.array(gamma_star)
    delta_star = np.array(delta_star)
    #for each batch
    for j, batch_idxs in enumerate(batch_info):

        dsq = np.sqrt(delta_star[j,:])
        dsq = dsq.reshape((len(dsq), 1))
        denom =  np.dot(dsq, np.ones((1, n_batches[j]))) #divide by sqrt delta_star
        numer = np.array(bayesdata[batch_idxs] - np.dot(batch_design.loc[batch_idxs], gamma_star).T) #subtract gamma_star

        bayesdata[batch_idxs] = numer / denom
    #multiply by square root of variance and add mean
    vpsq = np.sqrt(var_pooled).reshape((len(var_pooled), 1))
    bayesdata = bayesdata * np.dot(vpsq, np.ones((1, int(n_array)))) + stand_mean
    return bayesdata

def plot_histograms(df_norm,df,s_train,s_tissue):
    '''
    for each marker, return a histogram of trianing data and transformed data (df_norm)
    '''
    bins=50
    d_fig = {}
    for s_marker in df_norm.columns[df_norm.dtypes=='float64']:
        print(s_marker)
        fig,ax=plt.subplots(2,1,figsize = (3,4))
        for idxs, s_batch in enumerate(sorted(set(df_norm.batch))):
            df_batch = df_norm[(df_norm.batch==s_batch)].loc[:,s_marker] 
            if len(df_batch.dropna()) == 0:
                continue
            ax[0].hist(df.loc[df.index.str.contains(s_batch),s_marker],bins=bins,alpha=0.4, color=f'C{idxs}')
            ax[1].hist(df_batch,bins=bins,alpha=0.4, color=f'C{idxs}',label=s_batch)
            ax[0].set_yscale('log')
            ax[1].set_yscale('log')
            ax[0].set_title(f'{s_marker.split("_")[0]}: Raw Data')
            ax[1].set_title(f'{s_marker.split("_")[0]}: Combat')
            ax[1].legend()
        plt.tight_layout()
        plt.close()
        d_fig.update({s_marker:fig})
    return(d_fig)