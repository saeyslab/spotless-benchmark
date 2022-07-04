# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 14:53:36 2022

@author: Sai
"""

import pandas as pd
from PIL import Image, ImageDraw

#### INITIALIZATION #####
path = r"/home/chananchidas/spotless-benchmark/"

dataset_types = ["artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                  "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                  "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse"]
datasets = ['brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'pbmc', 'scc_p5']
ref=False
ref_text = "_ref" if ref else ""

# Each plot has a different layout :(
# x1, x2, y1, y2, col spacing, row spacing
#"prc": [111, 367, 88, 257, 16.5, 16],
dims_dict =  {"RMSE":[172, 555, 133, 400, 24.7, 23.5],
              "prc": [172, 555, 133, 400, 24.7, 23.5]}
#if ref:
#    dims_dict["precision"] = [111, 367, 87, 256, 16.5, 15]
    
colors_dict = {"cell2location":"#f8766d", "music":"#a3a500", "rctd":"#00bf7d",\
               "spotlight":"#00b0f6", "stereoscope":"#e76bf3", "tie":"#a1a1a1"}

##### MAIN CODE #####
for metric in dims_dict:
    # Read output from R
    df = pd.read_csv(path+ r'plots/faceted_plots/best_values/best_values_' + metric + '.tsv', sep="\t")
    
    # Create new section to count ties
    df['comb'] = df['dataset'] + "," + df['dataset_type']
    reps = df['comb'].value_counts()
    
    # Store best performing method (or tie) in matrix
    results = [[0]*8 for i in range(7)]
    for i, dataset in enumerate(datasets):
        for j, dt in enumerate(dataset_types):
            if reps[dataset+","+dt] > 1:
                results[i][j] = "tie"
            else:
                results[i][j] = df[(df['dataset']==dataset) & \
                                   (df['dataset_type']==dt)]["method"].values[0]
    
    # Calculate facet width and height
    d = dims_dict[metric]
    start_x, start_y = d[0], d[2]
    height, width = d[3]-d[2], d[1]-d[0]
    space_x, space_y = d[4], d[5]
    
    # Draw rectangles on image
    im = Image.open(path + r"plots/faceted_plots//" + metric + ref_text + '.png')
    im = im.convert("RGB")
    img1 = ImageDraw.Draw(im)  
    for i in range(7):
        for j in range(8):
            # Some math
            x = j*(width+space_x)+start_x
            y = i*(height+space_y)+start_y
            shape = [(x, y), (x+width, y+height)]
            
            # Outline is the best performing-method
            img1.rectangle(shape, outline=colors_dict[results[i][j]], width=7)
    im.save('/home/chananchidas/Pictures/' + metric + ref_text + '_boxed.png')
