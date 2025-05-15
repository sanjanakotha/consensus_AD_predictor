from tkinter import S
import warnings
from xmlrpc.client import boolean
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns
#sns.set_theme()
import protfasta
warnings.filterwarnings('ignore')
from matplotlib.ticker import PercentFormatter
from scipy.stats import ks_2samp
sns.set(style = "ticks", rc={'figure.figsize':(6.5,4)}, font_scale = 1.2)
import matplotlib.ticker as mtick
import math
from itertools import product
from matplotlib.colors import LogNorm

import AD_predictor_tools
import AD_comparison_tools

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

#Control colors:
LambertTFs_color = "Grey"
GSL_color = "Orange"
predictions_color = "Blue"
proteome_color = "Purple"

def add_AAcounts_and_allhydros(df, col_name = "ProteinWindowSeq", AAs=["A","I","L","M","V","F","W","Y","N","C","Q","S","T","D","E","R","H","K","G","P"]):
    for aa in AAs:
            df[aa]= df.apply(lambda row: row[col_name].count(aa),axis=1)
    df['AllHydros'] = df[["W","Y","F","L"]].sum(axis=1)

def add_charge(df):
    df["Charge"] = df["K"] + df["R"] - df["D"] - df["E"]

def add_allhydros(df):
    df['AllHydros'] = df["W"] + df["F"] + df["Y"] + df["L"]

def tile_and_plot(column, tiled_Lambert, tiled_GS_list, tiled_proteome = 0, figsize=(6,4), 
                  dpi=80, xtick_interval=10, xtick_labelsize=10, bins_interval=10, 
                  ylabel = True, path_to_save = False, title = True, 
                  legend = True, legend_outside = False, xlabel = False, alpha = 0.7, element = "bars", 
                  fill = True, first_list = "Lambert TFs", second_list = "Gold Standard List", 
                  third_list = "Proteome",
                  show = True, new_figure = True, use_matplotlib = False, mid_bin_label = True, custom_title = False):
    if new_figure:
        figure(figsize = figsize, dpi=dpi)

    if title:
        plt.title(column)
    
    if ylabel:
        plt.ylabel('Percent')

    if type(tiled_proteome) == int:
         bins=np.arange(min(min(tiled_Lambert[column]),
                                min(tiled_GS_list[column]))
                    ,max(max(tiled_Lambert[column]),
                                max(tiled_GS_list[column])) + bins_interval,
                    bins_interval)
    else:
        bins=np.arange(min(min(tiled_Lambert[column]),
                                min(tiled_GS_list[column]),
                                min(tiled_proteome[column]))
                    ,max(max(tiled_Lambert[column]),
                                max(tiled_GS_list[column]),
                                max(tiled_proteome[column])) + bins_interval,
                    bins_interval)

    def bins_labels(bins, **kwargs):
            bin_w = (max(bins) - min(bins)) / (len(bins) - 1)
            
            if type(tiled_proteome) == int:
                display_bins =np.arange(min(min(tiled_Lambert[column]),
                             min(tiled_GS_list[column]))
                   ,max(max(tiled_Lambert[column]),
                             max(tiled_GS_list[column])) + xtick_interval,
                   xtick_interval)
            else:
                display_bins = np.arange(min(min(tiled_Lambert[column]),
                                        min(tiled_GS_list[column]), 
                                        min(tiled_proteome[column]))
                                    ,max(max(tiled_Lambert[column]),
                                        max(tiled_GS_list[column]), 
                                        max(tiled_proteome[column])) + xtick_interval,
                                    xtick_interval)
            # print(len(display_bins))  
            # print(display_bins)
            # print(len(np.arange(min(bins)+bin_w/2, max(bins) + xtick_interval, xtick_interval)))
            if mid_bin_label:
                plt.xticks(np.arange(min(bins)+bin_w/2, max(bins) + xtick_interval, xtick_interval), display_bins, **kwargs)
            else:
                plt.xticks(np.arange(min(bins), max(bins) + xtick_interval, xtick_interval), display_bins, **kwargs)
            plt.xlim(bins[0], bins[-1])

    #tiled_Lambert[column].hist(bins=bins, alpha=alpha, grid=False, weights=np.ones(len(tiled_Lambert.index)) / len(tiled_Lambert.index),label="Lambert TFs", color = LambertTFs_color)
    #tiled_GS_list[column].hist(bins=bins, alpha=alpha, grid=False, weights=np.ones(len(tiled_GS_list.index)) / len(tiled_GS_list.index),label="Gold Standard List", color = GSL_color)
    if use_matplotlib:    
        matplotlib.rc_file_defaults()
        plt.hist(tiled_Lambert[column], bins=bins, alpha=alpha, density = True,label= first_list, color = LambertTFs_color)
        plt.hist(tiled_GS_list[column], bins=bins, alpha=alpha, density = True,label= second_list, color = GSL_color)

    else:
        sns.histplot(tiled_Lambert, element = element, fill = fill, x= column, bins = bins, weights=np.ones(len(tiled_Lambert.index)) / len(tiled_Lambert.index),label=first_list, color = LambertTFs_color, alpha = alpha)
        sns.histplot(tiled_GS_list, element = element, fill = fill, x= column, bins = bins, weights=np.ones(len(tiled_GS_list.index)) / len(tiled_GS_list.index),label=second_list, color = GSL_color, alpha = alpha)

    #tiled_Lambert["group"] = first_list
    #tiled_GS_list["group"] = second_list
    #plot_df = AD_comparison_tools.df_list_to_df([tiled_Lambert, tiled_GS_list])
    
    if type(tiled_proteome) != int:
        #tiled_proteome["group"] = "proteome"
        #lot_df =  AD_comparison_tools.df_list_to_df([plot_df, tiled_proteome])     
        if use_matplotlib:   
            plt.hist(tiled_proteome[column], bins=bins, alpha=alpha, density = True,label=third_list, color = proteome_color)

            sns.set_theme()
        else:
            sns.histplot(tiled_proteome, element = element, fill = fill, x= column, bins = bins, weights=np.ones(len(tiled_proteome.index)) / len(tiled_proteome.index),label="Proteome", color = proteome_color, alpha = alpha)

    #sns.histplot(data=plot_df, x=column, hue="species", multiple="stack")

    if xtick_labelsize:
        bins_labels(bins, fontsize = xtick_labelsize)
    else:
        bins_labels(bins)

    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals = 0))

    if xlabel:
        plt.xlabel(xlabel)
    
    if legend:
        if legend_outside: 
            plt.legend(loc=(1.04, 0.6))

        else:
            plt.legend()

    #plt.axvline(0, color = 'gray', linestyle = '--', lw = 1)

    if custom_title:
        plt.title(custom_title)
    
    if path_to_save:
        plt.savefig(path_to_save)

    if show:
        plt.show()



# Returns probability of each X, Y pair in tiled_df
def tiled_df_to_abundance_df(tiled_df, X_val = "Charge", Y_val = "AllHydros"):
    abundance_series = tiled_df.groupby([X_val, Y_val]).size()
    abundance_df = abundance_series.to_frame()
    abundance_df.reset_index(inplace=True)
    abundance_df = abundance_df.rename(columns = {0: "count"})
    abundance_df["log10count"] = [math.log10(abundance_df["count"][i]) for i in range(len(abundance_df.index))]
    return abundance_df

def add_Y_axis_val(tiled_df, AAs, y):
    Y_axis_val = tiled_df[AAs[0]].copy(deep = True)
    for AA in AAs[1:]:
        Y_axis_val += tiled_df[AA]
    tiled_df[y] = Y_axis_val

# Plots heatmap of tiled_df. 
    # label = 
    # c = log10count or count
def plot_one_heatmap(tiled_df, label, x = "Charge", c = "log10count", cmap = "Greys", AAs = ["W", "F", "Y", "L"], s =21, x_min = math.inf, x_max = -math.inf, y_min = math.inf, y_max = -math.inf, seaborn = True):
    add_Y_axis_val(tiled_df, AAs, "AllHydros")

    pivot_tbl =  pd.crosstab(tiled_df['Charge'], tiled_df['AllHydros'], dropna = False )
    pivot_tbl = pivot_tbl.T
    
    end_x_min = int(min(min(pivot_tbl.columns), x_min))
    end_x_max = int(max(max(pivot_tbl.columns), x_max))
    end_y_min = int(min(min(pivot_tbl.index), y_min))
    end_y_max = int(max(max(pivot_tbl.index), y_max))
    
    # Add columns to pivot table    
    for i in range(end_x_min, end_x_max + 1):
        if i not in pivot_tbl.columns:
            pivot_tbl[i] = np.zeros(len(pivot_tbl))
    
    #Add rows to pivot table
    for i in range(end_y_min, end_y_max + 1):
        if i not in pivot_tbl.index:
            pivot_tbl.loc[i] = np.zeros(len(pivot_tbl.columns))

    if seaborn: 
        pivot_tbl = pivot_tbl.sort_index()
        pivot_tbl = pivot_tbl.sort_index(axis = 1)
        
        ax = sns.heatmap(pivot_tbl, cmap = cmap, norm=LogNorm(), square = True)
        ax.invert_yaxis() 
        
        plt.xticks(rotation=0)
        plt.yticks(rotation=0)

        ax.set_xticks(ax.get_xticks()[::2])
        ax.set_yticks(ax.get_yticks()[::2])

        plt.xlabel('Charge')
        plt.ylabel('Number of ' + ','.join(AAs))
    
    else:
        abundance_df = tiled_df_to_abundance_df(tiled_df)
        cs1 = plt.scatter(data = abundance_df, x = x, y = "AllHydros", c = c, cmap = cmap,
                        alpha=1, label = label, marker='s', s=s, linewidths=0)
        plt.colorbar(cs1, shrink=1, aspect=20,label= c + '(' + label + ')',ticks=[0,1,2,3,4])
    
        if label:
            plt.legend()

    #return ax

def plot_new_heatmap(tiled_df, label, x = "Charge", y = "AllHydros", c = "log10count", cmap = "Greys"):
    plt.figure(figsize = (8, 3), dpi = 300)
    plot_one_heatmap(tiled_df, label, x, y, c, cmap)
    plt.show()

def plot_two_heatmaps(tiled_df1, tiled_df2, label1, label2, x = "Charge", y = "AllHydros", c = "log10count", cmap1 = "Greys", cmap2 = "Greys"):
    plt.figure(figsize = (8, 3), dpi = 300)
    plot_one_heatmap(tiled_df1, label1, x, y, c, cmap1)
    plot_one_heatmap(tiled_df2, label2, x, y, c, cmap2)
    plt.show()

def plot_preds(LowerCorner = [-13,7], UpperCorner = [-9,10], LowerCorner_slope1=0, 
            LowerCorner_slope2=0, UpperCorner_slope1=1, UpperCorner_slope2='inf', slope = 1, AAs = ["W", "F", "Y", "L"]):
 
    LambertDF = AD_predictor_tools.makeTilingDF("../data/LambertTFs.fasta")
    add_allhydros(LambertDF)
    
    LowerCorner_line1=(LowerCorner_slope1*(LambertDF.Charge-LowerCorner[0])-(LambertDF.AllHydros-LowerCorner[1]))<= 0

    LowerCorner_line2=LowerCorner_slope2*(LambertDF.Charge-LowerCorner[0])-(LambertDF.AllHydros-LowerCorner[1])<= 0

    if UpperCorner_slope1=="inf":
        UpperCorner_line1= (LambertDF.Charge<=UpperCorner[0])
    else:
        UpperCorner_line1=(UpperCorner_slope1*(LambertDF.Charge-UpperCorner[0])-(LambertDF.AllHydros-UpperCorner[1]))<= 0

    if UpperCorner_slope2=="inf":
        UpperCorner_line2= (LambertDF.Charge<=UpperCorner[0])
    else:
        UpperCorner_line2=UpperCorner_slope2*(LambertDF.Charge-UpperCorner[0])-(LambertDF.AllHydros-UpperCorner[1]) >= 0 


    BothSet= (LowerCorner_line1 ) & (LowerCorner_line2 ) & (UpperCorner_line1 ) & (UpperCorner_line2)

    plt.scatter(LambertDF[BothSet].Charge,LambertDF[BothSet].AllHydros,label='Above the line N = %s'%sum(BothSet))#, color = "#BADDE1")
    plt.scatter([LowerCorner[0],UpperCorner[0]],[LowerCorner[1],UpperCorner[1]],c='r')

    new_lower_corner_x=UpperCorner[0]-(UpperCorner[1]-LowerCorner[1])/slope
    new_lower_corner_y=LowerCorner[1]

    plt.ylabel('Number of '+','.join(AAs)), plt.xlabel('Charge'),plt.legend()

    plt.ylim([0,25])

# Draws line between two points
def plot_line(point1, point2, color = "b", lw = 3, alpha = 1, seaborn = False):
    x_values = [point1[0], point2[0]]
    y_values = [point1[1], point2[1]]
    
    plt.plot(x_values, y_values, color = color, lw = lw, alpha = alpha)

def return_line_y_coords(x_min, x_max, func):
        return [func(x) for x in np.arange(x_min, x_max)]

def plot_VP16_CITED2(c = "r", marker1 = "o", marker2 = "*", marker1size = 2.5, marker2size= 4, alpha = 1, seaborn = False):
    # plt.plot(-13, 7, c = c, marker = marker, markersize = markersize, alpha = alpha)
    # plt.plot(-9, 10, c = c, marker = marker, markersize = markersize, alpha = alpha)

    x1 = -12.5
    y1 = 7.5

    x2 = -8.5
    y2 = 10.5

    if seaborn:
        x1 += 39
        x2 += 39

    plt.plot(x1, y1, c = c, marker = marker1, markersize = marker1size, alpha = alpha)
    plt.plot(x2, y2, c = c, marker = marker2, markersize = marker2size, alpha = alpha)

def plot_boundary_lines(alpha = 1, lw = 3, seaborn = False, orig_predictor_only = False, with_EF = False, color = 'b', new_predictor = False):
    neg_thirteen = -13
    neg_twelve = -11
    fifteen = 20
    neg_nine = -8
    neg_twenty_five = -35
    ten = 10
    eleven = 11
    seven = 7

    if seaborn:
        neg_thirteen += 39
        neg_nine += 39
        neg_twenty_five += 39
        neg_twelve += 39

    point2 = [neg_nine, fifteen]
    point4 = [neg_nine, ten]
    point5 = [neg_twenty_five, seven]
    point7 = [neg_nine, seven]

        
    point1 = [neg_thirteen, fifteen]
    point3 = [neg_twenty_five,eleven]
    new_point4 = [neg_nine, eleven]
    point6 = [neg_thirteen, seven]

    lower_point7 = [neg_nine, seven - 1]
    right_point7 = [neg_nine + 1, seven - 1]
    right_point2 = [neg_nine + 1, fifteen]
    lower_point5 =  [neg_twenty_five, seven - 1]

    if new_predictor == False:
        

        

        plot_line(point5, [neg_twelve, seven], alpha = alpha, lw = lw, seaborn = seaborn, color = color)
        plot_line(point4, [neg_twelve, seven], alpha = alpha, lw = lw, seaborn = seaborn, color = color)
        plot_line(point4, point2, alpha = alpha, lw = lw, seaborn = seaborn, color = color)

        if  not orig_predictor_only:
            plot_line(point1, point6, alpha = alpha, lw = lw, seaborn = seaborn, color = color)
            plot_line(point2, point7, alpha = alpha, lw = lw, seaborn = seaborn, color = color)
            plot_line(point3, new_point4, alpha = alpha, lw = lw, seaborn = seaborn, color = color)
            plot_line(point5, point7, alpha = alpha, lw = lw, seaborn = seaborn, color = color)
        
        

        if with_EF:
            plot_line(lower_point7, point7, alpha = alpha, lw = lw, seaborn = seaborn, color = color)
            plot_line(right_point7, lower_point7, alpha = alpha, lw = lw, seaborn = seaborn, color = color)
            plot_line(right_point2, right_point7, alpha = alpha, lw = lw, seaborn = seaborn, color = color)
            plot_line(lower_point5, lower_point7, alpha = alpha, lw = lw, seaborn = seaborn, color = color)
    
    else:
        plot_line([neg_thirteen, 6], right_point7, alpha = alpha, lw = lw, seaborn = seaborn, color = color)
        plot_line([neg_thirteen, 6], point1, alpha = alpha, lw = lw, seaborn = seaborn, color = color)
        plot_line(right_point2, right_point7, alpha = alpha, lw = lw, seaborn = seaborn, color = color)



def plot_region_labels(seaborn = False, bbox={'facecolor': 'white', 'alpha': 0.8, 'pad': 0}, with_EF = False, fontsize = 12):
    if seaborn:
        correction = 40.5
    else:
        correction = 0
    plt.text(-12.5 + correction, 13, 'A', fontweight = "bold", c = "Black", bbox = bbox, fontsize = fontsize)
    plt.text(-12.5 + correction, 9, 'B', fontweight = "bold", c = "Black", bbox = bbox, fontsize = fontsize)
    plt.text(-17 + correction, 9, 'C', fontweight = "bold", c = "Black", bbox = bbox, fontsize = fontsize)
    plt.text(-10.6 + correction, 7.5, 'D', fontweight = "bold", c = "Black", bbox = bbox, fontsize = fontsize)

    if with_EF:
        plt.text(-17 + correction, 6.1, 'E', fontweight = "bold", c = "Black", bbox = bbox, fontsize = fontsize)
        plt.text(-9.4 + correction, 13, 'F', fontweight = "bold", c = "Black", bbox = bbox, fontsize = fontsize)

# Returns probability of each value in col_name in tiled_df
def tiled_df_to_abundance_df_one_col(tiled_df, col_name):
    abundance_series = tiled_df.groupby([col_name]).size()
    abundance_df = abundance_series.to_frame()
    abundance_df.reset_index(inplace=True)
    abundance_df = abundance_df.rename(columns = {0: "count"})
    abundance_df[col_name + "_prob"] = abundance_df["count"] / len(tiled_df.index)
    return abundance_df

def expected_abundance_df(tiled_proteome):
    Charge_df = tiled_df_to_abundance_df_one_col(tiled_proteome, "Charge")
    Charge_df = Charge_df.set_index("Charge")

    AllHydros_df = tiled_df_to_abundance_df_one_col(tiled_proteome, "AllHydros")
    AllHydros_df = AllHydros_df.set_index("AllHydros")

    charge_allhydros_pairs = list(product(Charge_df.index, AllHydros_df.index))
    charges = []
    allhydros = []
    expected_counts = []

    def expected_count(charge, allhydros):
        return Charge_df["Charge_prob"][charge] * AllHydros_df["AllHydros_prob"][allhydros] * len(tiled_proteome.index)
    
    for pair in charge_allhydros_pairs:
        charges.append(pair[0])
        allhydros.append(pair[1])
        expected_counts.append(expected_count(pair[0], pair[1]))

    expected_df = pd.DataFrame({"Charge": charges,
             "AllHydros": allhydros,
             "expected_count": expected_counts})
    expected_df["log10expected_count"] = np.log10(expected_df["expected_count"])
    return expected_df

def plot_marg_hists(tiled_df, label, color, x = "Charge", y = "AllHydros", c = "log10count", cmap = "Greys", expected_hist = False, rect = False):
    if (not expected_hist):
        abundance_df = tiled_df_to_abundance_df(tiled_df)
    else:
        abundance_df = expected_abundance_df(tiled_df)

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    if rect:
        fig = plt.figure(figsize=(9.8, 4.6), dpi = 300)
        spacing_x = 0.1 / 4.6
        spacing_y = 0.1 / 9.6

        rect_histx = [left, bottom + height + spacing_x, width, 1.5/4.6]
        rect_histy = [left + width + spacing_y, bottom, 1.5/9.6, height]
        rect_scatter = [left, bottom, width, height]

        s = 21 * (9.8 / 8) * (4.6 / 3)
    else:
        fig = plt.figure(figsize=(8, 8))

        rect_histx = [left, bottom + height + spacing, width, 0.2]
        rect_histy = [left + width + spacing, bottom, 0.2, height]
        rect_scatter = [left, bottom, width, height]
        
        s = 100

    

    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)

    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    cs1 = ax.scatter(data = abundance_df, x = x, y = y, c = c, cmap = cmap,
                    alpha=1,label= label,marker='s',s = s,linewidths=0)

    ax_histx.hist(tiled_df[x], bins= np.arange(min(tiled_df[x]), max(tiled_df[x])), color = color)
    ax_histy.hist(tiled_df[y], bins= np.arange(min(tiled_df[y]), max(tiled_df[y])), orientation='horizontal', color = color)

    ax.set_xlabel(x)
    ax.set_ylabel(y)

    plt.colorbar(cs1, shrink=1, aspect=20,label=label,ticks=[0,1,2,3,4],drawedges=False)

    plt.show()


def plot_two_marg_hists(tiled_Lambert, tiled_GS_list, label1, label2, cmap1 = "Greys", cmap2 = "Oranges", color1 = LambertTFs_color, color2= GSL_color, rect = False):
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    if rect:
        factor = 1.8
        box_width = 9.6
        x_dim = factor * box_width
        
        fig = plt.figure(figsize= (x_dim, 4.6), dpi = 300)
        spacing_x = 0.1 / 4.6
        spacing_y = 0.1 / x_dim

        width = width * (8 / x_dim)

        rect_histx = [left, bottom + height + spacing_x, width, 1.5/4.6]
        rect_histy = [left + width + spacing_y, bottom, 1.5/box_width, height]
        rect_scatter = [left, bottom, width, height]

        s = 21 * (box_width / 8) * (4.6 / 3)
    else:
        fig = plt.figure(figsize=(8, 8))

        rect_histx = [left, bottom + height + spacing, width, 0.2]
        rect_histy = [left + width + spacing, bottom, 0.2, height]
        rect_scatter = [left, bottom, width, height]
        
        s = 100
    
    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)

    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    min_x = min([min(tiled_Lambert["Charge"]), min(tiled_GS_list["Charge"])])
    max_x = max([max(tiled_Lambert["Charge"]), max(tiled_GS_list["Charge"])])

    min_y = min(min(tiled_GS_list["AllHydros"]), min(tiled_Lambert["AllHydros"]))
    max_y = max(max(tiled_GS_list["AllHydros"]), max(tiled_Lambert["AllHydros"]))
    
    LambertTFs_abundanceDF = tiled_df_to_abundance_df(tiled_Lambert)
    GSL_abundanceDF = tiled_df_to_abundance_df(tiled_GS_list)

    cs1 = ax.scatter(data = LambertTFs_abundanceDF, x = "Charge", y = "AllHydros", c = "log10count", cmap = cmap1,
                    alpha=1,label=label1,marker='s',s=s,linewidths=0)
    cs2 = ax.scatter(data = GSL_abundanceDF, x = "Charge", y = "AllHydros", c = "log10count", cmap = cmap2,
                    alpha=1,label=label2,marker='s',s=s,linewidths=0)

    ax_histx.hist(tiled_Lambert["Charge"], bins= np.arange(min_x, max_x), color = color1, alpha = 0.5, density = True)
    ax_histx.hist(tiled_GS_list["Charge"], bins= np.arange(min_x, max_x), color = color2, alpha = 0.5, density = True)

    ax_histy.hist(tiled_Lambert["AllHydros"], bins= np.arange(min_y, max_y), orientation='horizontal', color = color1, alpha = 0.5, density = True)
    ax_histy.hist(tiled_GS_list["AllHydros"], bins= np.arange(min_y, max_y), orientation='horizontal', color = color2, alpha = 0.5, density = True)

    ax.set_xlabel("Charge")
    ax.set_ylabel("AllHydros")

    plt.colorbar(cs1, shrink=1, aspect=20,label=label1,ticks=[0,1,2,3,4],drawedges=False, pad = 0.4/factor)
    plt.colorbar(cs2, shrink=1, aspect=20,label=label2,ticks=[0,1,2,3,4],drawedges=False, pad = 0.1/factor)

    plt.show()
