import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve, roc_curve, roc_auc_score, auc
import numpy as np

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D

sns.set_theme(rc={'figure.figsize':(6,4), 'figure.dpi' : 300})
sns.set_style("white")

lambert_TFs = pd.read_csv("../data/LambertTFs.csv", index_col = 0)

all_known_ADs = pd.read_csv("../output/known_ADs_considering_isoforms.csv")
# Only keep ADs on Lambert TFs
lambert_TFs["uniprotID"] = lambert_TFs["GeneName"].str.split("|").str[1]
known_ADs = all_known_ADs[all_known_ADs["TileType"] == "TF"]
known_ADs = known_ADs[known_ADs["uniprotID"].isin(lambert_TFs["uniprotID"])]
known_ADs.head(3)


def add_overlap_status(df, domain_list, col_type, min_overlap = 32, col_name_suffix =  "_suffic_overlap"):
    # Keep rows of adhunter with annots
    # Calculate interval overlap
    # If there is overlap of at least 80% (so >= 32 residues) then there is sufficient overlap
    df_with_annots = pd.merge(df, domain_list, on = "uniprotID")
    df_with_annots["overlap_length"] =  df_with_annots[['End', 'annot_End']].min(axis=1) - df_with_annots[['Start', 'annot_Start']].max(axis=1)
    df_with_annots[col_type+ "_suffic_overlap"] = df_with_annots["overlap_length"] >= min_overlap
    #display(df_with_annots)
    
    # Keep first of each tile, such that if there is an overlappign tile it is saved
    df_with_annots = df_with_annots.sort_values(by = col_type+ "_suffic_overlap", ascending = False)
    df_with_annots = df_with_annots.drop_duplicates(subset = ["tile"], keep = "first")
    #display(df_with_annots)
    
    # Add info back to full tile table (includes TFs with no annotated ADs)
    #display(df)
    df = pd.merge(df, df_with_annots[["tile", col_type+ "_suffic_overlap"]], how = "left")
    df = df.fillna(False)
    return df

def process_overlaps(df):
    # Keep rows with at least one annotations
    return_df = df[(df["AD_suffic_overlap"]) | (df["RD_suffic_overlap"])]
    
    # Keep rows where AD and RD aren't both true
    return_df = return_df[~(return_df["AD_suffic_overlap"] & return_df["RD_suffic_overlap"])]

    return_df["active"] = (return_df["AD_suffic_overlap"]) & (~return_df["RD_suffic_overlap"])
    return return_df

def process_overlaps_with_dbd(df):
    # Keep rows with at least one annotations
    return_df = df[(df["AD_suffic_overlap"]) | (df["RD_suffic_overlap"]) | (df["DBD_suffic_overlap"])]
    
    # Keep rows where only one annot type
    return_df["annot_count"] = return_df["AD_suffic_overlap"].astype(int) + return_df["RD_suffic_overlap"].astype(int) + return_df["DBD_suffic_overlap"].astype(int)
    return_df = return_df[return_df["annot_count"] == 1]

    return_df["active"] = (return_df["AD_suffic_overlap"]) & (~return_df["RD_suffic_overlap"]) & (~return_df["DBD_suffic_overlap"])
    return return_df

def print_overlap_counts(pred_overlap, pred_ad_rd, pred_ad_rd_dbd):
    print("Total predictions: " + str(len(pred_overlap)))
    print("AD or RD overlapping: " + str(len(pred_ad_rd)))
    print("AD or RD or DBD overlapping: " + str(len(pred_ad_rd_dbd)))

# Helper functions
def formatting():
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.gca().set_aspect('equal', adjustable='box')
    sns.despine()

def plot_prc(df, pred_name, color="b", active_col_name="active"):
    precision, recall, thresholds = precision_recall_curve(df[active_col_name], df[pred_name])
    plt.plot(recall, precision, color=color)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    formatting()
    return auc(recall, precision)

def plot_roc(df, pred_name, color="b", first=True, text=True, active_col_name="active"):
    fpr, tpr, thresholds = roc_curve(df[active_col_name], df[pred_name])
    plt.plot(fpr, tpr, color=color)
    plt.xlabel("FPR")
    plt.ylabel("TPR")

    # Calculate AUROC and display it on the plot
    auroc = roc_auc_score(df[active_col_name], df[pred_name])
    if text:
        y = 0.15 if first else 0.05
        plt.text(x=0.95, y=y, s=f"AUROC = {auroc:.2f}", fontsize='small', va="bottom", ha="right", color=color)

    # Find the best threshold (closest to the top-left corner)
    best_idx = np.argmax(tpr - fpr)
    best_threshold = thresholds[best_idx]
    print(f"Best threshold for {pred_name}: {best_threshold:.4f}")

    formatting()
    return auroc

def add_custom_legend(color_dict, ax, bbox_to_anchor=(1, 1)):
    custom_lines = [Line2D([0], [0], markersize=2, color=c, lw=4) for c in color_dict.values()]
    ax.legend(custom_lines, color_dict.keys(), bbox_to_anchor=bbox_to_anchor, frameon=False, fontsize='x-small', handlelength=0.5)

def plot_prc_roc(df, pred_name, active_col_name='active'):
    activator_color = sns.color_palette('colorblind')[1]
    all_color = sns.color_palette('colorblind')[0]
    
    # Create two subplots sharing the y-axis
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(7, 4))
    activator_TF = df[df["uniprotID"].isin(known_ADs["uniprotID"])]

    # Plot Precision-Recall on the first subplot
    plt.sca(ax1)  # Switch to ax1
    plt.plot([0, 1], [0, 1], 'k--', color="gray", lw=1)
    plot_roc(df, pred_name, all_color, first=True, active_col_name=active_col_name)
    plot_roc(activator_TF, pred_name, activator_color, first=False, active_col_name=active_col_name)

    # Plot ROC curve on the second subplot
    plt.sca(ax2)  # Switch to ax2
    all_TF_random_luck = sum(df[active_col_name]) / len(df)
    plt.axhline(all_TF_random_luck, linestyle='--', color=all_color, lw=1)

    activ_TF_random_luck = sum(activator_TF[active_col_name]) / len(activator_TF)
    plt.axhline(activ_TF_random_luck, linestyle='--', color=activator_color, lw=1)

    plot_prc(df, pred_name, all_color, active_col_name=active_col_name)
    plot_prc(activator_TF, pred_name, activator_color, active_col_name=active_col_name)

    add_custom_legend({"All TFs": all_color, "Activator TFs": activator_color}, ax=ax2)
    plt.subplots_adjust(wspace=0.2, hspace=0)


def plot(pred_output, col_name, min_overlap = 32):
    
    # Adding uniprotID to adhunter to merge with known ADs
    pred_output["uniprotID"] = pred_output["GeneName"].str.split("|").str[1]
    pred_output = pred_output.rename(columns = {"StartPosition" : "Start", "EndPosition" : "End", "ProteinWindowSeq" : "tile"})
    pred_output

    # Preparing known AD coords for merge
    known_AD_coords = known_ADs[["uniprotID", "Start", "End"]]
    known_AD_coords = known_AD_coords.rename(columns = {"Start" : "annot_Start", "End" : "annot_End"})
    known_AD_coords

    pred_output_overlap = add_overlap_status(pred_output, known_AD_coords, "AD", min_overlap=min_overlap)
    pred_output_overlap

    known_RDs = pd.read_csv("../../SFARI/data/tycko_soto_delrosso_RD_coordinate_data_updated_12-2-24.csv", header = None)
    known_RDs[1] = known_RDs[1].str.strip("[").str.strip("]").str.split(")")
    known_RDs = known_RDs.explode(1)
    known_RDs = known_RDs[known_RDs[1] != ""]
    known_RDs["annot_Start"] = known_RDs[1].str.extract(r'\((\d+),').astype(int)
    known_RDs["annot_End"] = known_RDs[1].str.extract(r', (\d+)$').astype(int)
    known_RDs["uniprotID"] = known_RDs[0]
    known_RDs = known_RDs[["uniprotID", "annot_Start", "annot_End"]]
    known_RDs

    pred_output_overlap = add_overlap_status(pred_output_overlap, known_RDs, "RD", min_overlap=min_overlap)
    pred_output_overlap

    known_DBDs = pd.read_csv("../../SFARI/output/lambert_TFs_10-21-24_with_DBD_coords.csv", index_col = 0)
    known_DBDs = known_DBDs.dropna()
    known_DBDs["DBD_coords_merged"] = known_DBDs["DBD_coords_merged"].str.split("], ")
    known_DBDs = known_DBDs.explode("DBD_coords_merged")
    known_DBDs["annot_Start"] = known_DBDs["DBD_coords_merged"].str.extract(r'\[(\d+),').astype(int)
    known_DBDs["annot_End"] = known_DBDs["DBD_coords_merged"].str.extract(r', (\d+)').astype(int)
    # known_RDs["uniprotID"] = known_RDs[0]
    # known_RDs = known_RDs[["uniprotID", "annot_Start", "annot_End"]]
    known_DBDs["uniprotID"] = known_DBDs["id"].str.split("|").str[1]
    known_DBDs

    pred_output_overlap = add_overlap_status(pred_output_overlap, known_DBDs, "DBD", min_overlap = min_overlap)
    pred_output_overlap

    processed_output = process_overlaps_with_dbd(pred_output_overlap)

    sns.set_style('ticks')
    sns.set_context('talk')
    plot_prc_roc(processed_output, col_name)
