import os
import skimage
from skimage import io, segmentation, morphology, measure
import numpy as np
import tifffile
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
#import bioformats 
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
