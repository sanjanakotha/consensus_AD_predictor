import pandas as pd
import numpy as np

#Returns dataframe with variations with a higher proportion of predictions overapping with GS than orig_GSoverlap
def higher_gs_prop_preds(outputfilepath="../../PredictionADs_ToShare/Output/Varying save_listope/output_100_perms_interval_0.5_LambertTFs_withaLine.csv",
                                    filepath_before_variation="../../PredictionADs_ToShare/Output/Varying save_listope/",
                                    filepath_after_variation="output_predicted_activation_domains_LambertTFs_withaLine.csv",
                                    orig_GSoverlap=26/144,
                                    variation_column_name="Thresholdsave_listope"):   

    output_df=pd.read_csv(outputfilepath).drop(axis=1,labels="Unnamed: 0")
    output_df=output_df[output_df["ProportionPredictedAD_GoldStandardOverlaps"]>=orig_GSoverlap]
    return_df=[]

    for variation in output_df[variation_column_name]:
        variation_df=pd.read_csv(filepath_before_variation+str(variation)+filepath_after_variation,index_col=0)
        variation_df["Description"]=variation_column_name+str(variation)+" with Higher GS Prop"
        return_df.append(variation_df)

    print("There are "+str(len(return_df))+" variations with a higher proportion of entries on the gold standard list.")

    n_preds=0
    for df in return_df:
        n_preds+=len(df.index)
        
    print("There are "+str(n_preds)+" sequences to test (before checking for overlaps).")

    return df_list_to_df(return_df)


#Take in a list of dataframes and appends all of them into one dataframe which it returns
def df_list_to_df(df_list, note_list=None, note_list_col_name=None):
    if type(df_list)==list:
        return_df=pd.DataFrame()
        for i in range(len(df_list)):
            if note_list:
                df_list[i][note_list_col_name]=note_list[i]
            return_df=return_df.append(df_list[i], ignore_index=True)
        return return_df
    else:
        print("input is not a list of dataframes!")
        return df_list


#Takes in one column of a dataframe and returns a list with the number of identical entries in the column for each entry
def count_occurrences_in_series(df_column):
    return_list=[]
    df_column_list=df_column.tolist()
    for entry in df_column:
        return_list.append(df_column_list.count(entry))
    return return_list  

#Adds a column to a dataframe with the number of times a value in a column appeared
def add_occurrence_count(df,column_name):
    df["times_"+column_name+"_predicted"]=count_occurrences_in_series(df[column_name])

# Drops duplicate ProteinRegionSeqs in the input dataframe; updates Start, End, Length; and merges the original descriptions.
# def remove_duplicate_sequences(df):

# Helper function, source: https://stackoverflow.com/questions/27182137/check-if-two-lines-each-with-start-end-positions-overlap-in-python
def are_separate(r, s):
  (a, b) = r
  (x, y) = s
  return b < x or y < a

# Function to check if a dataframe contains a prediction, given a uniprotID and a start/end
# returns True or False
def contains_prediction(pred_df_row, compare_to_df, ID_col_name="uniprotID", start_col_name="Start", end_col_name="End"):

    #rows of compare_to_df with same uniprotID
    compare_to_df=compare_to_df[compare_to_df[ID_col_name]==pred_df_row[ID_col_name]]
    
    #is there a row where compare to pred? 
    for start_col_val, end_col_val in zip(compare_to_df[start_col_name], compare_to_df[end_col_name]):
        if not are_separate((start_col_val,end_col_val),(pred_df_row[start_col_name],pred_df_row[end_col_name])):
            return True
    
    return False

#Adds a column to pred_df with whether each row is contained in compare_to_df or not.
def add_col_contains_prediction(pred_df, compare_to_df, ID_col_name="uniprotID", start_col_name="Start", end_col_name="End", result_col_name="contained_in_df2"):
    results=[]
    for i in range(len(pred_df.index)):
        results.append(contains_prediction(pred_df.iloc[i], compare_to_df))
    pred_df[result_col_name]=results

#returns a dataframe with a status for each prediction indicating if it is in the set compareTo. 
#If "PredictionCounts">0, the prediction is in compareTo set
#compareTo must have uniprotIDs in the GeneName.
def overlapstatus(predictionDF,compareTo):
    
    #Testing if predictionDF has the necessary information.
    if set(['Start','End','uniprotID']).issubset(predictionDF.columns):
        print("Columns of DataFrame 1 look good.")
    elif "GeneName" in predictionDF.columns:
        if "|" in predictionDF["GeneName"][0]:
            predictionDF["uniprotID"]=predictionDF.apply(lambda row: row['GeneName'].split('|')[1],axis=1)
            print("Columns of DataFrame 1 look good after getting uniprotIDs from GeneName.")
    else:
        print("There is not enough information in the the first dataframe to run this function.")
        print("The first dataframe needs a Start and End column and either a uniprotID column or a GeneName from which you can get a uniprotID.")

    #Testing if compareTo has the necessary information.
    if set(['Start','End','uniprotID']).issubset(compareTo.columns):
        print("Columns of DataFrame 2 look good.")
    elif "|" in compareTo["GeneName"][0]:
        compareTo["uniprotID"]=compareTo.apply(lambda row: row['GeneName'].split('|')[1],axis=1)
        print("Columns of DataFrame 2 look good after getting uniprotID from GeneName.")
    else:
        print("There is not enough information in the the second dataframe to run this function.")
        print("The second dataframe needs a Start and End column and either a uniprotID column or a GeneName from which you can get a uniprotID.")
    
    #---------------------------------------------------------
    
    Nregions = len(predictionDF.index)

    if "Activity_mean" in compareTo.columns:
        Activity_mean_list= np.zeros(Nregions)
        OverlapsKRAB_list=np.zeros(Nregions)
    
    ## actual predictions
    predictionCounts,overlapcounter = np.zeros(Nregions),0
        
    for i,entry_uniprotID,entry_Start,entry_End in zip(np.arange(0,len(predictionDF.index)),predictionDF["uniprotID"],predictionDF["Start"],predictionDF["End"]): # note I choose the random regions in the same order each time, but the TFs are randomized so it should be fine    ## if TF has a known AD
        ## note, this only generates a random region of the TF if the selected TF has a known AD. Otherwise, the chance of overlap is Zero and to save time,I do not choose a random region
        if entry_uniprotID in compareTo.uniprotID.values:
            predictedStart = entry_Start
            predictedEnd = entry_End
            
            # assess if sampled regions overlap know ADs
            indx = compareTo.uniprotID==entry_uniprotID
            TFwithKnownADs =  compareTo[indx]
            
            if "Activity_mean" in compareTo.columns:
                # Because a TF can have multiple AD entries, cycle through these entries
                for KnownStart,KnownEnd,Activity_mean,OverlapsKRAB in zip(TFwithKnownADs["Start"],TFwithKnownADs["End"],TFwithKnownADs["Activity_mean"],TFwithKnownADs["OverlapsKRAB"]):
        #                 print(entry_GeneName)
        #                 print('%i - %i Known'%(KnownStart,KnownEnd))
        #                 print('%i - %i Tested'%(predictedStart,predictedEnd)) 
                    if (predictedStart<=KnownEnd)&(predictedEnd>KnownStart):
                        overlapcounter +=1
                        predictionCounts[i]+=1
                        Activity_mean_list[i]=Activity_mean
                        OverlapsKRAB_list[i]=OverlapsKRAB
            else:
                # Because a TF can have multiple AD entries, cycle through these entries
                for KnownStart,KnownEnd in zip(TFwithKnownADs["Start"],TFwithKnownADs["End"]):
        #                 print(entry_GeneName)
        #                 print('%i - %i Known'%(KnownStart,KnownEnd))
        #                 print('%i - %i Tested'%(predictedStart,predictedEnd)) 
                    if (predictedStart<=KnownEnd)&(predictedEnd>KnownStart):
                        overlapcounter +=1
                        predictionCounts[i]+=1

                
    Npredictions =np.sum(predictionCounts>0)
    #print('There are %i predictions total. %i were tested. %i were not tested.'%(Nregions,Npredictions,Nregions-Npredictions))
    
    overlap_status_df=predictionDF
    overlap_status_df["OverlapStatus"]=predictionCounts
    if "Activity_mean" in compareTo.columns:
        overlap_status_df["Activity_mean"]=Activity_mean_list
        overlap_status_df["OverlapsKRAB"]=OverlapsKRAB_list

    #returns [predictions that are in the tested set,  
    #number of predictions not in the set]
    
    #print("There are "+str(len(predictionDF.index))+" predictions. "+str(len(tested_predictions_df))+" are in the set with uniprotIDs.")
    
    return overlap_status_df

# both columns need dataframes with labels uniprotID, Start, and End
def compare_two_predictors(our_preds_df, other_preds_df, comp_list = pd.read_csv("../data/Gold Standard AD List.csv"), pred1_name = "Our Predictor", pred2_name = "Other Predictor", comp_list_name = "GSL"):
    our_preds = our_preds_df.copy(deep=True)
    other_preds = other_preds_df.copy(deep=True)
    
    # Number of our_preds that overlap with the GSL
    add_col_contains_prediction(our_preds, comp_list, result_col_name = "GSL")
    our_preds_overlap_count = sum(our_preds["GSL"])
    our_preds_unique_GSL_overlap_count = unique_GSL_count(our_preds, GSL = comp_list)
    
    # Number of other_preds that overlap with the GSL
    add_col_contains_prediction(other_preds, comp_list, result_col_name = "GSL")
    other_preds_overlap_count = sum(other_preds["GSL"])
    other_preds_unique_GSL_overlap_count = unique_GSL_count(other_preds, GSL = comp_list)

    # Our predictions that the other predictor also makes
    add_col_contains_prediction(our_preds, other_preds, result_col_name = "other_preds")
    our_preds_and_other = our_preds[our_preds["other_preds"]]
    # Number of our_preds confirmed by other_preds that overlap with the GSL
    add_col_contains_prediction(our_preds_and_other, comp_list, result_col_name = "GSL")
    our_preds_and_other_overlap_count = sum(our_preds_and_other["GSL"])
    our_preds_and_other_unique_GSL_overlap_count = unique_GSL_count(our_preds_and_other, GSL = comp_list)
    
    # The other predictor's predictions that our predictor also makes
    add_col_contains_prediction(other_preds, our_preds, result_col_name = "our_preds")
    other_preds_and_us = other_preds[other_preds["our_preds"]]
    # Number of other_preds confirmed by our_preds that overlap with the GSL
    add_col_contains_prediction(other_preds_and_us, comp_list, result_col_name = "GSL")
    other_preds_and_us_overlap_count = sum(other_preds_and_us["GSL"])
    other_preds_and_us_unique_GSL_overlap_count = unique_GSL_count(other_preds_and_us, GSL = comp_list)
    
    pred1_confirmed_by_pred2_name = pred1_name + " confirmed by " + pred2_name
    pred2_confirmed_by_pred1_name = pred2_name + " confirmed by " + pred1_name
    
    predicted_by = [pred1_name, pred2_name, pred1_confirmed_by_pred2_name, pred2_confirmed_by_pred1_name]
    number_preds = [len(our_preds.index), len(other_preds.index), len(our_preds_and_other.index), len(other_preds_and_us.index)]
    GSL_overlap_count = [our_preds_overlap_count, other_preds_overlap_count, our_preds_and_other_overlap_count, other_preds_and_us_overlap_count]
    unique_GSL_overlap_count = [our_preds_unique_GSL_overlap_count, other_preds_unique_GSL_overlap_count, our_preds_and_other_unique_GSL_overlap_count, other_preds_and_us_unique_GSL_overlap_count]
    
    return_df = pd.DataFrame({"predicted_by": predicted_by,
                        "number_preds": number_preds,
                        "num_preds_on_" + comp_list_name: GSL_overlap_count,
                        "num_"+comp_list_name+"_entries_in_preds": unique_GSL_overlap_count})

    return_df["prop_preds_on_" + comp_list_name] = return_df["num_preds_on_" + comp_list_name] / return_df["number_preds"]

    return_df["prop_"+comp_list_name+"_entries_in_preds"] = return_df["num_"+comp_list_name+"_entries_in_preds"] / return_df["number_preds"]
    
    our_preds = our_preds.drop(columns = ["GSL", "other_preds"])
    our_preds_and_other = our_preds_and_other.drop(columns = ["other_preds"])
    other_preds = other_preds.drop(columns = ["GSL", "our_preds"])
    other_preds_and_us = other_preds_and_us.drop(columns = ["our_preds"])

    return return_df

def unique_GSL_count(preds, GSL = pd.read_csv("../data/Gold Standard AD List.csv")):
    if type(preds) == str:
        predictionDF = pd.read_csv(preds)
    else:
        predictionDF = preds
    

    if "uniprotID" not in predictionDF.columns:
        predictionDF["uniprotID"] = predictionDF.apply(lambda row: row['GeneName'].split('|')[1],axis=1)
    #GSL = pd.read_csv(GSL_filepath)
    add_col_contains_prediction(GSL, predictionDF, result_col_name = "containedInPreds")
    return_value = sum(GSL["containedInPreds"])
    GSL.drop(["containedInPreds"], axis = 1)

    return return_value

def GSL_count(orig_preds, GSL_filepath = "../data/newGSL.csv"):
    # Adding a column to original predictions with status of whether it overlaps with GSL
    GSL = pd.read_csv(GSL_filepath)
    add_col_contains_prediction(orig_preds, 
                                GSL, 
                                ID_col_name="uniprotID", 
                                start_col_name="Start", 
                                end_col_name="End", 
                                result_col_name="contained_in_GSL")
    ours_GSL = len(orig_preds[orig_preds["contained_in_GSL"]==True].index)
    ours_len = len(orig_preds.index)
    print(str(round(ours_GSL/ours_len, 3))+", or, "+str(ours_GSL)+
      " out of "+str(ours_len)+" predictions made by us are on the gold standard list.")
    return ours_GSL


# +
# To merge lists of ADs (ex. the GSL and Soto)

# Helper function: returns merged row
def return_merged_row(uniprotID, df):
    # Only look at rows with the same uniprot ID
    same_uniprotID_rows = df[df["uniprotID"] == uniprotID]
    same_uniprotID_rows = same_uniprotID_rows.sort_values(by = "Start")
    
    # Final dataframe columns
    new_starts = []
    new_ends = []
    genes = []
    AD_names = []
    references = []
    orig_lists = []
    
    # Current row's values
    curr_start = -1
    curr_end = -1
    curr_genes = []
    curr_AD_names = []
    curr_references = []
    curr_lists = []
    
    for i in same_uniprotID_rows.index:
        # Merge current row with next row
        if curr_end >= same_uniprotID_rows.loc[i]["Start"]:
            curr_end = max(curr_end, same_uniprotID_rows.loc[i]["End"])
            curr_genes.append(same_uniprotID_rows.loc[i]["GeneName"])
            curr_AD_names.append(same_uniprotID_rows.loc[i]["AD name"])
            curr_references.append(same_uniprotID_rows.loc[i]["Reference"])
            curr_lists.append(same_uniprotID_rows.loc[i]["List"])
        
        # Don't merge current row with next row
        else: 
            new_starts.append(curr_start)
            new_ends.append(curr_end)
            genes.append(" / ".join(set([c.strip() for c in curr_genes])))
            
            curr_AD_names = [str(c) for c in curr_AD_names]
            AD_names.append(" / ".join(curr_AD_names))
            
            curr_references = [str(c) for c in curr_references]
            references.append(" / ".join(curr_references))
            
            curr_lists = [str(c) for c in curr_lists]
            orig_lists.append(" / ".join(curr_lists))
            
            curr_start = same_uniprotID_rows.loc[i]["Start"]
            curr_end = same_uniprotID_rows.loc[i]["End"]
            
            curr_genes = [same_uniprotID_rows.loc[i]["GeneName"]]
            curr_AD_names = [same_uniprotID_rows.loc[i]["AD name"]]
            curr_references = [same_uniprotID_rows.loc[i]["Reference"]]
            curr_lists = [same_uniprotID_rows.loc[i]["List"]]
    
    # Append the last values
    new_starts.append(curr_start)
    new_ends.append(curr_end)
    
    genes.append(" / ".join(set([c.strip() for c in curr_genes])))
    
    curr_AD_names = [str(c) for c in curr_AD_names]
    AD_names.append(" / ".join(curr_AD_names))
    
    curr_references = [str(c) for c in curr_references]
    references.append(" / ".join(curr_references))
    
    curr_lists = [str(c) for c in curr_lists]
    orig_lists.append(" / ".join(curr_lists))
    
    # Remove the first (because it is just -1 or "")
    new_starts = new_starts[1:]
    new_ends = new_ends[1:]
    genes = genes[1:]
    AD_names = AD_names[1:]
    references = references[1:]
    orig_lists = orig_lists[1:]
    
    return pd.DataFrame({"GeneName": genes,
                         "AD name": AD_names,
                         "Start": new_starts,
                        "End": new_ends,
                        "uniprotID": uniprotID,
                         "Reference": references,
                         "List": orig_lists
                        })

# Input: list of dataframes with "GeneName", "AD name", "Start", "End", "uniprotID", "Reference" columns
# Output: one dataframe with merged entries (merges entries with any overlap)
def return_merged_list(df_list):
    both_lists = pd.concat(df_list)
    both_lists = both_lists.reset_index(drop = True)
    dfs = []
    for uniprotID in both_lists["uniprotID"].unique():
        dfs.append(return_merged_row(uniprotID, both_lists))
    new_GSL = pd.concat(dfs)
    return new_GSL.sort_values(by = "uniprotID").reset_index(drop = True)
# -


