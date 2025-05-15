import pandas as pd
import numpy as np

# This function takes in a dataframe with TFs split into multiple rows, then merges them back into one row. 
# based on function in notebook: Running ADPred on Lambert TFs
def merge_tiled_sequences(split_LambertTFs):
    
    merged_LambertTFs=pd.DataFrame({"GeneName":[],"ProteinSeq":[],"Length":[], "combined position_wise_prob_adpred":[]})

    i=0

    while i < len(split_LambertTFs.index):

        #If the next GeneName is not identical or there is no next row, add the row as is. 
        if i == (len(split_LambertTFs.index) - 1) or split_LambertTFs["GeneName"][i] != split_LambertTFs["GeneName"][i+1]: 
            arr = activity_string_to_array(split_LambertTFs["position_wise_prob_adpred"][i])
            merged_LambertTFs.loc[len(merged_LambertTFs.index)] = [split_LambertTFs["GeneName"][i],
                                                                   split_LambertTFs["ProteinWindowSeq"][i],
                                                                   len(split_LambertTFs["ProteinWindowSeq"][i]),
                                                                   str(arr)]


        #Else, if the next GeneName is identical, 
        else:
            prot_seq=split_LambertTFs["ProteinWindowSeq"][i][:-15]
            position_wise_prob_adpred = activity_string_to_array(split_LambertTFs["position_wise_prob_adpred"][i])[:-15]

            #Keep concatenating the protein seqs while the next GeneName is the same. update i each time. 
            while i != (len(split_LambertTFs.index) - 1) and split_LambertTFs["GeneName"][i] == split_LambertTFs["GeneName"][i+1]:
                prot_seq+=split_LambertTFs["ProteinWindowSeq"][i+1][15:][:-15]  
                position_wise_prob_adpred += activity_string_to_array(split_LambertTFs["position_wise_prob_adpred"][i+1])[15:][:-15]
                i+=1

            prot_seq+=split_LambertTFs["ProteinWindowSeq"][i][-15:]
            position_wise_prob_adpred += activity_string_to_array(split_LambertTFs["position_wise_prob_adpred"][i])[-15:]

            merged_LambertTFs.loc[len(merged_LambertTFs.index)] = [split_LambertTFs["GeneName"][i-1], 
                                                                 prot_seq, 
                                                                 len(prot_seq),
                                                                 position_wise_prob_adpred] 

        i+=1
        
    return merged_LambertTFs


def activity_string_to_array(prob_string):
    a = prob_string.replace('\n','')
    a = a.replace('[','')
    a = a.replace(']','').split(',')
    b = []
    for x in a:
        if x != "":
            b.append(float(x))
    return b

# arr = integer array position_wise_prob_adpred
def return_pred_df(arr, min_prob, min_length, ProteinSeq, GeneName):
    prob_arr = [a >= min_prob for a in arr]

    start_indices = []
    end_indices = []

    i = 0
    
    while i < (len(arr) - 1):
        start = i
        temp = [arr[i]]
        add = prob_arr[i]
        while(i < len(arr) - 1 and prob_arr[i] and prob_arr[i+1]):
            temp.append(arr[i+1])
            i += 1
            add = True
        if(len(temp) >= min_length and add):
            start_indices.append(start)
            end_indices.append(i)
        i += 1

    if (len(arr) != len(ProteinSeq)):
        print("ERROR: arr length != sequence length")
        print("Array length: "+str(len(arr)))
        print("ProteinSeq length: "+str(len(ProteinSeq)))
        print("GeneName")
    protein_seqs = [ProteinSeq[s:e+1] for s, e in zip(start_indices, end_indices)]
    
    return_df = pd.DataFrame({"GeneName": GeneName,
                            "ProteinRegionSeq": protein_seqs,
                            "Start": start_indices,
                            "End": end_indices})
        
    return return_df

