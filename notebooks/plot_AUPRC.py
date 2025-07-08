import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve, roc_curve, roc_auc_score, auc
import numpy as np

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
import AD_comparison_tools

sns.set_theme(rc={'figure.figsize':(6,4), 'figure.dpi' : 300})
sns.set_style("white")

lambert_TFs = pd.read_csv("../../SFARI/output/lambert_TFs_10-21-24_with_DBD_coords.csv", index_col = 0)

# Only keep ADs on Lambert TFs
all_known_ADs = pd.read_csv("../../SFARI/output/known_ADs_considering_isoforms_and_canonical_with_alerasool.csv")
lambert_TFs["uniprotID"] = lambert_TFs["id"].str.split("|").str[1]
known_ADs = all_known_ADs[all_known_ADs["TileType"] == "TF"]
known_ADs = known_ADs[known_ADs["uniprotID"].isin(lambert_TFs["uniprotID"])]
known_ADs.head(3)

# Preparing known AD coords for merge
known_AD_coords = known_ADs[["uniprotID", "Start", "End"]]
known_AD_coords = known_AD_coords.rename(columns = {"Start" : "annot_Start", "End" : "annot_End"})

# Preparing known RD coords for merge
known_RDs = pd.read_csv("../../SFARI/data/tycko_soto_delrosso_RD_coordinate_data_updated_12-2-24.csv", header = None)
known_RDs[1] = known_RDs[1].str.strip("[").str.strip("]").str.split(")")
known_RDs = known_RDs.explode(1)
known_RDs = known_RDs[known_RDs[1] != ""]
known_RDs["annot_Start"] = known_RDs[1].str.extract(r'\((\d+),').astype(int)
known_RDs["annot_End"] = known_RDs[1].str.extract(r', (\d+)$').astype(int)
known_RDs["uniprotID"] = known_RDs[0]
known_RDs = known_RDs[["uniprotID", "annot_Start", "annot_End"]]

# Preparing known DBD coords for merge
known_DBDs = pd.read_csv("../../SFARI/output/lambert_TFs_10-21-24_with_DBD_coords.csv", index_col = 0)
known_DBDs = known_DBDs.dropna()
known_DBDs["DBD_coords_merged"] = known_DBDs["DBD_coords_merged"].str.split("], ")
known_DBDs = known_DBDs.explode("DBD_coords_merged")
known_DBDs["annot_Start"] = known_DBDs["DBD_coords_merged"].str.extract(r'\[(\d+),').astype(int)
known_DBDs["annot_End"] = known_DBDs["DBD_coords_merged"].str.extract(r', (\d+)').astype(int)
known_DBDs["uniprotID"] = known_DBDs["id"].str.split("|").str[1]


def add_overlap_status(df, domain_list, col_type, min_overlap, col_name_suffix =  "_suffic_overlap"):
    # Keep rows of adhunter with annots
    # Calculate interval overlap
    # If there is overlap of at least min_overlap then there is sufficient overlap
    df_with_annots = pd.merge(df, domain_list, on = "uniprotID")
    display(df_with_annots)
    df_with_annots["overlap_length"] =  df_with_annots[['End', 'annot_End']].min(axis=1) - df_with_annots[['Start', 'annot_Start']].max(axis=1)
    df_with_annots[col_type+ "_suffic_overlap"] = df_with_annots["overlap_length"] >= min_overlap
    #display(df_with_annots)
    
    # Keep first of each tile, such that if there is an overlappign tile it is saved
    df_with_annots = df_with_annots.sort_values(by = col_type+ "_suffic_overlap", ascending = False)
    df_with_annots = df_with_annots.drop_duplicates(subset = ["tile"], keep = "first")
    #display(df_with_annots)
    
    # Add info back to full tile table (includes TFs with no annotated ADs)
    #ns.histplot(df_with_annots[df_with_annots["overlap_length"] > 0]["overlap_length"])
    df = pd.merge(df, df_with_annots[["tile", col_type+ "_suffic_overlap"]], how = "left")
    df = df.fillna(False)
    return df

# def process_overlaps(df):
#     # Keep rows with at least one annotations
#     return_df = df[(df["AD_suffic_overlap"]) | (df["RD_suffic_overlap"])]
    
#     # Keep rows where AD and RD aren't both true
#     return_df = return_df[~(return_df["AD_suffic_overlap"] & return_df["RD_suffic_overlap"])]

#     return_df["active"] = (return_df["AD_suffic_overlap"]) & (~return_df["RD_suffic_overlap"])
#     return return_df

def process_overlaps_with_ad_and_dbd(df):
    # Keep rows with at least one annotation
    return_df = df[(df["AD_suffic_overlap"]) | (df["DBD_suffic_overlap"])]
    
    # Keep rows where only one annot type
    return_df["annot_count"] = return_df["AD_suffic_overlap"].astype(int) + return_df["DBD_suffic_overlap"].astype(int)
    return_df = return_df[return_df["annot_count"] == 1]

    return_df["active"] = (return_df["AD_suffic_overlap"]) & (~return_df["DBD_suffic_overlap"])
    return return_df

def process_overlaps_with_rd_and_dbd(df):
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
def formatting(ax=None):
    if ax is None:
        ax = plt.gca()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal', adjustable='box')
    sns.despine(ax=ax)

def plot_prc(df, pred_name, color="b", active_col_name="active", ax=None, subtract_prevalence=True):
    precision, recall, thresholds = precision_recall_curve(df[active_col_name], 
                                                           df[pred_name])
    if ax is None:
        ax = plt.gca()

    ax.plot(recall, precision, color=color)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    formatting(ax)

    if subtract_prevalence:
        return auc(recall, precision) - sum(df[active_col_name]) / len(df)
    else:
        return auc(recall, precision)

def plot_roc(df, pred_name, color="b", first=True, text=True, active_col_name="active", ax=None, custom_threshold=None, return_threshold=False):
    fpr, tpr, thresholds_roc = roc_curve(df[active_col_name], df[pred_name])
    if ax is None:
        ax = plt.gca()
    ax.plot(fpr, tpr, color=color)
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")

    # Calculate AUROC and display it on the plot
    auroc = roc_auc_score(df[active_col_name], df[pred_name])
    if text:
        y = 0.15 if first else 0.05
        ax.text(x=0.95, y=y, s=f"AUROC = {auroc:.2f}", fontsize='small', va="bottom", ha="right", color=color)

    # Find the best threshold based on TPR - FPR
    best_idx = np.argmax(tpr - fpr)
    best_threshold = thresholds_roc[best_idx]
    best_fpr = fpr[best_idx]
    best_tpr = tpr[best_idx]

    # Plot the best threshold point on ROC
    ax.plot(best_fpr, best_tpr, 'o', color=color, markersize=10)

    #     # If custom threshold is given, compute and plot it
    if custom_threshold is not None:
        tp = len(df[df["active"] & (df[pred_name] > custom_threshold)])
        fp = len(df[~df["active"] & (df[pred_name] > custom_threshold)])
        tn = len(df[~df["active"] & ~(df[pred_name] > custom_threshold)])
        fn = len(df[df["active"] & ~(df[pred_name] > custom_threshold)])
        
        # # Example usage: print the values (or handle them as needed)
        # print(f"Custom Threshold: {custom_threshold}, TP: {tp}, FP: {fp}, TN: {tn}, FN: {fn}")

        # Compute TPR and FPR
        custom_tpr = tp / (tp + fn)
        custom_fpr = fp / (fp + tn)

        print(custom_tpr, custom_fpr)
        print()

    #     # Plot custom threshold point
    #     ax.plot(custom_fpr, custom_tpr, 's', color='black', markersize=8, label=f"Custom Threshold ({custom_threshold:.2f})")



    # print(f"Best threshold for {pred_name} (F1 score): {best_threshold:.4f}")
    # print(f"Best F1 score for {pred_name}: {f1_scores[best_idx]:.4f}")
    formatting(ax)
    if return_threshold:
        return auroc, best_threshold
    else:
        print("Best threshold for " + pred_name + ":" + str(best_threshold))
        return auroc
    
def add_custom_legend(color_dict, ax, bbox_to_anchor=(1, 1)):
    custom_lines = [Line2D([0], [0], markersize=2, color=c, lw=4) for c in color_dict.values()]
    ax.legend(custom_lines, color_dict.keys(), bbox_to_anchor=bbox_to_anchor, frameon=False, fontsize='x-small', handlelength=0.5)

def plot_prc_roc_for_pair(df, pred_name, active_col_name='active', ax=None):
    activator_color = sns.color_palette('colorblind')[1]
    all_color = sns.color_palette('colorblind')[0]
    
    if ax is None:
        fig, ax = plt.subplots(1, 2, sharey=True, figsize=(7, 4))

    ax1, ax2 = ax

    activator_TF = df[df["uniprotID"].isin(known_ADs["uniprotID"])]

    # Plot ROC curve on the first subplot
    ax1.plot([0, 1], [0, 1], 'k--', color="gray", lw=1)
    plot_roc(df, pred_name, all_color, first=True, active_col_name=active_col_name, ax=ax1)
    plot_roc(activator_TF, pred_name, activator_color, first=False, active_col_name=active_col_name, ax=ax1)

    # Plot Precision-Recall on the second subplot
    all_TF_random_luck = sum(df[active_col_name]) / len(df)
    ax2.axhline(all_TF_random_luck, linestyle='--', color=all_color, lw=1)

    activ_TF_random_luck = sum(activator_TF[active_col_name]) / len(activator_TF)
    ax2.axhline(activ_TF_random_luck, linestyle='--', color=activator_color, lw=1)

    plot_prc(df, pred_name, all_color, active_col_name=active_col_name, ax=ax2)
    plot_prc(activator_TF, pred_name, activator_color, active_col_name=active_col_name, ax=ax2)

    add_custom_legend({"All TFs": all_color, "Activator TFs": activator_color}, ax=ax2)
    plt.subplots_adjust(wspace=0.2, hspace=0)

# def return_processed_output(tiled_TFs, min_overlap):
#     pred_output_overlap = add_overlap_status(tiled_TFs, known_AD_coords, "AD", min_overlap=min_overlap)
#     pred_output_overlap = add_overlap_status(pred_output_overlap, known_DBDs, "DBD", min_overlap=min_overlap)
#     processed_output = process_overlaps_with_ad_and_dbd(pred_output_overlap)
#     return processed_output

def plot_one_predictor_pair(pred_output, col_name, min_overlap):

    # Adding uniprotID to adhunter to merge with known ADs
    if "uniprotID" not in pred_output.columns:
        pred_output["uniprotID"] = pred_output["GeneName"].str.split("|").str[1]
    pred_output = pred_output.rename(columns = {"StartPosition" : "Start", "EndPosition" : "End", "ProteinWindowSeq" : "tile"})

    pred_output_overlap = add_overlap_status(pred_output, known_AD_coords, "AD", min_overlap=min_overlap)
    #pred_output_overlap = add_overlap_status(pred_output_overlap, known_RDs, "RD", min_overlap=min_overlap)
    pred_output_overlap = add_overlap_status(pred_output_overlap, known_DBDs, "DBD", min_overlap = min_overlap)
    #display(pred_output_overlap)

    processed_output = process_overlaps_with_ad_and_dbd(pred_output_overlap)
    #display(processed_output)

    sns.set_style('ticks')
    sns.set_context('talk')
    plot_prc_roc_for_pair(processed_output, col_name)

    return processed_output

def return_processed_output(pred_output, min_overlap):
    # Adding uniprotID to adhunter to merge with known ADs
    if "uniprotID" not in pred_output.columns:
        pred_output["uniprotID"] = pred_output["GeneName"].str.split("|").str[1]
    pred_output = pred_output.rename(columns = {"StartPosition" : "Start", "EndPosition" : "End", "ProteinWindowSeq" : "tile"})

    pred_output_overlap = add_overlap_status(pred_output, known_AD_coords, "AD", min_overlap=min_overlap)
    pred_output_overlap = add_overlap_status(pred_output_overlap, known_DBDs, "DBD", min_overlap = min_overlap)

    processed_output = process_overlaps_with_ad_and_dbd(pred_output_overlap)
    
    return processed_output



def return_tp_fp_v2(pred_name, binarized_summary_table):
    tp_count = len(binarized_summary_table[binarized_summary_table["active"] & binarized_summary_table[pred_name]])
    fp_count = len(binarized_summary_table[~binarized_summary_table["active"] & binarized_summary_table[pred_name]])
    return tp_count, fp_count

def eval_model2(model_name, binarized_summary_table):
    #preds = model_dict[model_name]
    TF_lim_preds = binarized_summary_table[binarized_summary_table["uniprotID"].isin(known_ADs["uniprotID"])]
    print(len(TF_lim_preds))
    num_preds = sum(TF_lim_preds[model_name])
    positive_benchmark_count = sum(TF_lim_preds["active"])
    negative_benchmark_count = sum(~TF_lim_preds["active"])
    TP, FP = return_tp_fp_v2(model_name, TF_lim_preds)
    return model_name, num_preds, positive_benchmark_count, negative_benchmark_count, TP, FP


def plot_random_expec(df, color, ax):
    all_TF_random_luck = sum(df["active"]) / len(df)
    ax.axhline(all_TF_random_luck, linestyle = '--', color = color, lw = 1)
    return all_TF_random_luck
