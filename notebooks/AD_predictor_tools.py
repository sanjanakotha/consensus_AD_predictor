from turtle import color
import numpy as np
import protfasta
import pandas as pd
import matplotlib.pyplot as plt
import random
from matplotlib.pyplot import figure
import PlottingTools
import seaborn as sns
sns.set_theme()
sns.set(style = "ticks", rc={'figure.figsize':(6.5,4)}, font_scale = 1.2)


# extract 39 AA windows from all TFs in a FASTA file
def makeTilingDF(inputfilename,window_size = 39, window_spacing = 1,AAs = ['W','F','Y','M','L','Q']):
    import os.path
    from os import path
    import re

    tilingDFs = os.listdir("../data/TilingDFs/")
    df=inputfilename.split("/")[-1].split(".")[0]
    #regexp = re.compile(df+".+size_"+str(window_size)+".+space_"+str(window_spacing)+r'.+AAs.+'+'.+'.join(AAs)+'.+')

    AA_string=""
    for AA in AAs:
        AA_string+="(?=.*"+AA+")"

    regexp = re.compile(df + ".+size_" + str(window_size) + ".+space_" + str(window_spacing) + r'.+AAs.+' + AA_string + '.+')

    newlist = list(filter(regexp.match, tilingDFs)) 

    filepath = "../data/TilingDFs/" + inputfilename.split("/")[-1].split(".")[0] + "_size_" + str(window_size) + "_space_" + str(window_spacing) + "_AAs_" + ",".join(AAs) + ".csv"
    
    if bool(newlist):
        print("Using existing Tiling DF at ../data/TilingDFs/" + newlist[0])
        return pd.read_csv("../data/TilingDFs/" + newlist[0], index_col = 0, low_memory = True)
    else:       
        print("Creating new Tiling DF at ",filepath)     
        ProteomeDict=protfasta.read_fasta(inputfilename, return_list=False, 
                                    invalid_sequence_action='ignore',
                                    duplicate_record_action='ignore',
                                    duplicate_sequence_action='ignore',
                                expect_unique_header=False)
        
        TilingWindows,ProteinNames,IUpredScores,Npos,Cpos = [],[],[],[],[]

        for gene in ProteomeDict:
                # take tilingwindows
                line = ProteomeDict[gene]
                if len(line)<window_size:
                    subset = line
                    Npos.append(0),Cpos.append(len(line))
                    TilingWindows.append(subset)
                    ProteinNames.append(gene)
                else:
                    for i in np.arange(0,len(line)-window_size+window_spacing,window_spacing):
                        subset = line[i:i+window_size]
                        Npos.append(i),Cpos.append(i+len(subset))
                        TilingWindows.append(subset)
                        ProteinNames.append(gene)

        print("Window Size = %i  and Window spacing = %i" % (window_size, window_spacing))
        print("Number of Tiling Windows: %i" %len(TilingWindows))

        ProteomeDF = pd.DataFrame({'ProteinWindowSeq':TilingWindows,'GeneName':ProteinNames,
                                'StartPosition':Npos,'EndPosition':Cpos})
        
        AAs = AAs + ['K', 'R', 'D', 'E']
        #count important AAs
        for aa in AAs:
            # ProteomeDF[aa]= ProteomeDF.apply(lambda row: row['ProteinWindowSeq'].count(aa),axis=1)
            ProteomeDF[aa] = ProteomeDF["ProteinWindowSeq"].str.count(aa)
        
        #compute net charge
        # ProteomeDF['Charge'] = ProteomeDF.apply(lambda row: row['ProteinWindowSeq'].count('K')+row['ProteinWindowSeq'].count('R')-row['ProteinWindowSeq'].count('D')-row['ProteinWindowSeq'].count('E'),axis=1)
        ProteomeDF['Charge'] = ProteomeDF['K'] +  ProteomeDF['R'] - ProteomeDF['D'] - ProteomeDF['E']

        ProteomeDF.to_csv(filepath)

        return ProteomeDF

def makeTilingDF_fromDF(inputDF, window_size = 39, window_spacing = 1, AAs = ['W','F','Y','M','L','Q'], col_name = "Sequence", add_gene_id = True):
    import random
    import os
    
    if add_gene_id:
        inputDF["GeneID"]=[">sp|"+str(a)+"|"+str(b)+"|"+str(c)+"|"+str(d)
                    for a,b,c,d in zip(inputDF["uniprotID"],
                                            inputDF["GeneName"],
                                            inputDF["Start"],
                                            inputDF["End"])]
    else:
        inputDF["GeneID"] = inputDF["uniprotID"]
        
    inputDF_dict=pd.Series(inputDF[col_name].values,index=inputDF.GeneID).to_dict()

    inputfilename = '../data/temp'+ str(random.getrandbits(128)) + '.fasta'

    protfasta.write_fasta(inputDF_dict, inputfilename)

    rV = makeTilingDF(inputfilename, window_size, window_spacing, AAs)

    filepath = "../data/TilingDFs/" + inputfilename.split("/")[-1].split(".")[0] + "_size_" + str(window_size) + "_space_" + str(window_spacing) + "_AAs_" + ",".join(AAs) + ".csv"
    os.remove(filepath)

    return rV

def df_to_fasta(inputDF, filename, seq_col_name):
    inputDF["GeneID"]=["sp|"+str(a)+"|"+str(b)+"|"+str(c)+"|"+str(d)
                   for a,b,c,d in zip(inputDF["uniprotID"],
                                        inputDF["GeneName"],
                                        inputDF["Start"],
                                        inputDF["End"])]

    inputDF_dict=pd.Series(inputDF[seq_col_name].values,index=inputDF.GeneID).to_dict()
    protfasta.write_fasta(inputDF_dict, filename)

# also make a DF with full length protein sequences
def makeFullLengthProteinDF(inputfilename):


    ProteinSeqs = protfasta.read_fasta(inputfilename, return_list=True, 
                                 invalid_sequence_action='ignore',
                                  duplicate_record_action='ignore',
                                  duplicate_sequence_action='ignore',
                                expect_unique_header=False)
    ProteomeDF=pd.DataFrame(ProteinSeqs).rename(columns={0:"GeneName",1:"AAseq"})
    print('There are %i proteins'%len(ProteomeDF.index)) 
    return ProteomeDF

def ThresholdProteome(inputfilename,
                        ProteomeDF,
                        LowerCorner,
                        UpperCorner,
                        LowerCorner_slope1, LowerCorner_slope2,
                        UpperCorner_slope1, UpperCorner_slope2,
                        AAs,
                        slope,
                        window_size,
                        propset =['Charge','AllHydros'],
                        plot_new_red_points = True):

    fig = plt.figure(figsize=(8,3),dpi=300)
    AllHydros = ProteomeDF[AAs].sum(axis=1)
    ProteomeDF['AllHydros'] = AllHydros
    new_heatmap = False

    # if(inputfilename == "../data/LambertTFs.fasta" and set(AAs) == set(["W", "F", "Y", "L"])):
    #     if (propset[0]=='Charge')&(propset[1]=='AllHydros'):
    #         overlayPredictionsOverAllTiles_v2(propset)
    # else: 
    tiled_df = makeTilingDF(inputfilename, window_size = window_size)
    ax = PlottingTools.plot_one_heatmap(AAs = AAs, tiled_df = tiled_df, label = 'All Human TF regions')

    if LowerCorner_slope1 == "inf":
        LowerCorner_line1 = ProteomeDF.Charge >= LowerCorner[0]
    else:
        LowerCorner_line1=(LowerCorner_slope1*(ProteomeDF.Charge-LowerCorner[0])-(ProteomeDF.AllHydros-LowerCorner[1]))<= 0

    LowerCorner_line2=LowerCorner_slope2*(ProteomeDF.Charge-LowerCorner[0])-(ProteomeDF.AllHydros-LowerCorner[1])<= 0
    
    if UpperCorner_slope1=="inf":
        UpperCorner_line1= (ProteomeDF.Charge<=UpperCorner[0])
    else:
        UpperCorner_line1=(UpperCorner_slope1*(ProteomeDF.Charge-UpperCorner[0])-(ProteomeDF.AllHydros-UpperCorner[1]))<= 0

    if UpperCorner_slope2=="inf":
        UpperCorner_line2= (ProteomeDF.Charge<=UpperCorner[0])
    else:
        UpperCorner_line2 = UpperCorner_slope2 * (ProteomeDF.Charge-UpperCorner[0]) - (ProteomeDF.AllHydros-UpperCorner[1]) >= 0 

    #LineValue = slope*(ProteomeDF.Charge-UpperCorner[0])-(ProteomeDF.AllHydros-UpperCorner[1])
    
    BothSet= (LowerCorner_line1 ) & (LowerCorner_line2 ) & (UpperCorner_line1 ) & (UpperCorner_line2)

    #BothSet = LowerCorner_line2 & (ProteomeDF.Charge<=UpperCorner[0]) & (ProteomeDF.AllHydros>=LowerCorner[1])
    #BothSet = (LineValue<=0) & (ProteomeDF.Charge<=UpperCorner[0]) & (ProteomeDF.AllHydros>=LowerCorner[1])
    #BothSet = (ProteomeDF.AllHydros<=UpperCorner[1]) & (ProteomeDF.AllHydros>=LowerCorner[1])
    #BothSet = (ProteomeDF.Charge<=UpperCorner[0]) & (ProteomeDF.Charge>=LowerCorner[0])


    seaborn_adjustment = min(ProteomeDF.Charge)
    plt.scatter(ProteomeDF[BothSet].Charge - seaborn_adjustment + 0.5, ProteomeDF[BothSet].AllHydros + 0.5,label='Above the line N = %s'%sum(BothSet))#, color = "#BADDE1")
    
    if plot_new_red_points:
        plt.scatter([LowerCorner[0] - seaborn_adjustment+ 0.5, UpperCorner[0] - seaborn_adjustment+ 0.5], [LowerCorner[1] + 0.5, UpperCorner[1] +0.5],color='r')

    print(ProteomeDF[BothSet].AllHydros)

    print("seaborn adjustment: " + str(seaborn_adjustment))

    new_lower_corner_x=UpperCorner[0]-(UpperCorner[1]-LowerCorner[1])/slope
    new_lower_corner_y=LowerCorner[1]

    #plt.scatter(new_lower_corner_x,new_lower_corner_y,color='orange')

    #Plotting lines:
    
    # #Lower corner line 1:
    # m=LowerCorner_slope1
    # x=ProteomeDF.Charge
    # plt.plot(x,m*(x-LowerCorner[0])+LowerCorner[1],color='#6DD3CE')

    # #Lower corner line 2:
    # m=LowerCorner_slope2
    # plt.plot(x,m*(x-LowerCorner[0])+LowerCorner[1],color='#C8E9A0')

    # #Upper corner line 1:
    
    # if UpperCorner_slope2=="inf":
    #     plt.axvline(x=UpperCorner[0],color='#F7A278')
    # else:
    #     m=UpperCorner_slope1
    #     plt.plot(x,m*(x-UpperCorner[0])+UpperCorner[1],color='#F7A278')

    # #Lower corner line 2:
    # if UpperCorner_slope2=="inf":
    #     plt.axvline(x=UpperCorner[0],color='#A13D63')
    # else:
    #     m=UpperCorner_slope2
    #     plt.plot(x,m*(x-UpperCorner[0])+UpperCorner[1],color='#A13D63')

    
    
    plt.ylabel('Number of '+','.join(AAs)), plt.xlabel('Charge'),plt.legend()

    plt.ylim([0,25])


    return BothSet

# ------------------------------------------------------------------------------------------

def MakeMesh(tempDF,propset,auto_propspan=True,min_prop1=-39,max_prop1=21,min_prop2=0,max_prop2=20):
    if auto_propspan:
        Prop1Span = np.arange(min(tempDF[propset[0]]),max(tempDF[propset[0]]),1)
        Prop2Span = np.arange(min(tempDF[propset[1]]),max(tempDF[propset[1]]),1)
    else:
        Prop1Span = np.arange(min_prop1,max_prop1,1)# hard coded to match shape of precomputed heatmap
        Prop2Span = np.arange(min_prop2,max_prop2,1) # hard coded to match shape of precomputed heatmap
    PredictedCounts = np.zeros((len(Prop1Span),len(Prop2Span)))

    Prop1SpanMesh, Prop2SpanMesh = np.meshgrid(Prop1Span,Prop2Span)

    # calculate number of tiles at each point on the grid.
    for i, p1 in enumerate(Prop1Span):
        for j,  p2 in enumerate(Prop2Span):
            indx = (tempDF[propset[0]]==p1)&(tempDF[propset[1]]==p2)
            PredictedCounts[i,j] = sum(indx)
    return Prop1SpanMesh, Prop2SpanMesh, PredictedCounts

def overlayPredictionsOverAllTiles_v2(propset =['Charge','AllHydros'],flip=False):
    Prop1Span = np.arange(-39,21,1)# hard coded to match shape of precomputed heatmap
    Prop2Span = np.arange(0,20,1) # hard coded to match shape of precomputed heatmap
    Prop1SpanMesh, Prop2SpanMesh = np.meshgrid(Prop1Span,Prop2Span)
    ## Human TFs 
    Counts_df=pd.read_csv("../data/HumanTFTiles_Counts_v2.csv")
    Counts=Counts_df.values
    Counts=Counts[:,1:]
    Z =np.log10(Counts.transpose())

    if flip:
        Y,X= Prop1SpanMesh, Prop2SpanMesh,
    else:
        X,Y= Prop1SpanMesh, Prop2SpanMesh,

    cs1  = plt.scatter(X,Y, c=Z, cmap="Greys", alpha=1,label='All Human TF regions',marker='s',s=21,linewidths=0)

    # Add a color bar which maps values to colors.
    plt.colorbar(cs1, shrink=1, aspect=20,label='log10(Abundance)',ticks=[0,1,2,3,4],drawedges=False)

# ------------------------------------------------------------------------------------------

# This is the helper function that runs prediction function. 
def maskproteome(inputfilename,AAs,Sequences_to_Test, LowerCorner_slope1, LowerCorner_slope2,
                        UpperCorner_slope1, UpperCorner_slope2, slope, window_size, lineparameters=[[-13,7],[-9,10]],propset=["Charge","AllHydros"], plot_new_red_points = True):
    #TEMPindx = ThresholdProteome(thresholdType,Sequences_to_Test,lineparameters[0],lineparameters[1],lineparameters[2],AAs)
    TEMPindx = ThresholdProteome(inputfilename,Sequences_to_Test,lineparameters[0],lineparameters[1], LowerCorner_slope1, LowerCorner_slope2,
                        UpperCorner_slope1, UpperCorner_slope2, AAs, slope, window_size, propset=["Charge","AllHydros"], plot_new_red_points=plot_new_red_points)

    MaskedProteomeDF = Sequences_to_Test[TEMPindx]
    if (len(MaskedProteomeDF) > 0):
        print('There are %i regions of length %i AA as extreme or more than this AD'%(len(MaskedProteomeDF),len(MaskedProteomeDF.ProteinWindowSeq.values[0]))) 
    else:
        print('There are %i regions as extreme or more than this AD'%(len(MaskedProteomeDF)))
    
    print('These regions come from %i proteins' % len(set(MaskedProteomeDF.GeneName)))

    ## aggregate the regions
    for gene in set(MaskedProteomeDF.GeneName):
        indx = MaskedProteomeDF.GeneName == gene
        tempDF = MaskedProteomeDF[indx]
        start, end = min(tempDF.StartPosition),max(tempDF.EndPosition)
        #print('%s \t%i AA \tfrom %i - %i'%(gene,end-start+1,start,end)) 
    return TEMPindx
# ------------------------------------------------------------------------------------------

#returns a file with format: "output_predictedADs_inputfilename_s:slope_h:allhydros_c:charge_composition_tl:tilinglength_ws:window_spacing"
def return_exportfilename(folder_name="predictions/",
                        inputfilename="'../data/LambertTFs.fasta", 
                        slope=1,
                        lower_corner_c=-13,lower_corner_h=7,
                        upper_corner_c=-9,upper_corner_h=10,
                        LowerCorner_slope1=1, LowerCorner_slope2=1,
                        UpperCorner_slope1=1, UpperCorner_slope2=1,
                        composition=["W","F","Y","L"],
                        window_size=39,
                        window_spacing=1,
                        propset=["Charge","AllHydros"]):
    return (folder_name+
            inputfilename.rpartition(".")[0].rpartition("/")[-1].zfill(10)
            +"_s_"+str(slope).zfill(3) #done
            +"_lcc_"+str(lower_corner_c).zfill(3)
            +"_lch_"+str(lower_corner_h).zfill(3)
            +"_ucc_"+str(upper_corner_c).zfill(3)
            +"_uch_"+str(upper_corner_h).zfill(3)
            +"_lcs1_"+str(LowerCorner_slope1).zfill(3)
            +"_lcs2_"+str(LowerCorner_slope2).zfill(3)
            +"_lcs1_"+str(UpperCorner_slope1).zfill(3)
            +"_ucs2_"+str(UpperCorner_slope2).zfill(3)
            +"_comp_"+"".join(composition) #done
            +"_tl_"+str(window_size).zfill(3) #done 
            +"_ws_"+str(window_spacing).zfill(3)
            +"_ps1_"+str(propset[0])
            +"_ps2_"+str(propset[1]))
            #+"_charge_only") #done

# Function to combine the Predicted Tiles into predicted ADs--look for overlapping tiles
def AggregateTilesIntoPredictedADs(inputfilename, Sequences_to_Test_DF,TEMPindx,exportfilename):

    CombinedPredictions = Sequences_to_Test_DF[TEMPindx]

    genenames, starts,ends, seqs,regionlengths,RegionType = [],[],[],[],[],[]
    for gene in set(CombinedPredictions.GeneName):
        indx = CombinedPredictions.GeneName == gene
        tempDF = CombinedPredictions[indx]
        tempDF = tempDF.reset_index()
        start, end = min(tempDF.StartPosition)+1,max(tempDF.EndPosition)
        # test if there are missing tiles or multiple regions
        if (((max(tempDF.StartPosition)-min(tempDF.StartPosition))-(sum(indx)))>39): #check if tiles overlap
            # tiles DO NOT overlap, use hand called regions
            previousstart,previousend = tempDF.StartPosition[0],tempDF.EndPosition[0]
            for i, StartPosition,EndPosition in zip(range(len(tempDF.index)),tempDF["StartPosition"],tempDF["EndPosition"]):
                if StartPosition < previousend+39:
                    previousend = EndPosition
                else:
                    starts.append(previousstart)
                    ends.append(previousend)
                    regionlengths.append(previousend-previousstart)
                    genenames.append(gene)
                    RegionType.append('Prediction')
                    previousstart,previousend = StartPosition,EndPosition
            # add last predicted region
            starts.append(previousstart)
            ends.append(EndPosition)
            regionlengths.append(EndPosition-previousstart)
            genenames.append(gene)
            RegionType.append('Prediction')
            if 0==1: # debugging
                print(gene +' had more than one region')
                print(tempDF[['GeneName','StartPosition','EndPosition']])
        else:
            starts.append(start)
            ends.append(end)
            regionlengths.append(end-start)
            genenames.append(gene)
            RegionType.append('Prediction')
    print("\n---\nThere are %i predicted candidate AD regions on %i TFs" %(len(starts),len(set(genenames))))
    # make a DF with geneNames and start and End positions of the overlapping windows
    CandidateADsToTest = pd.DataFrame({'GeneName':genenames,'Start':starts,'End':ends,'Length':regionlengths,'RegionType':RegionType})

    #add a column for full length AA seq
    FullLengthProteinsTested = makeFullLengthProteinDF(inputfilename)

    #print(FullLengthProteinsTested)
    tempDict = dict(zip(FullLengthProteinsTested.GeneName,FullLengthProteinsTested.AAseq))
    #print(len(tempDict))
    #print(tempDict.keys())
    tempSeries = CandidateADsToTest.GeneName.astype(str)
    #print(tempSeries)
    CandidateADsToTest['FullProteinSeq'] = tempSeries.map(tempDict)
    #print(CandidateADsToTest['FullProteinSeq'])

    ## pull out regions of these TFs
    # CandidateADsToTest['Region2Test'] = CandidateADsToTest.apply(lambda row: row.FullProteinSeq[row.Start:row.End], axis=1)
    #correct for counting starting at 0 in python and 1 in real world
    ProteinRegionSeqs =[]
    for i, Start, End, FullProteinSeq in zip(range(len(CandidateADsToTest.index)),CandidateADsToTest["Start"],CandidateADsToTest["End"],CandidateADsToTest["FullProteinSeq"]):
        start =Start
        end = End
        #display(CandidateADsToTest)
        #print("seq")
        #print(FullProteinSeq)
        Region = FullProteinSeq[start:end]
        ProteinRegionSeqs.append(Region)
    CandidateADsToTest['ProteinRegionSeq']=ProteinRegionSeqs
    CandidateADsToTest = CandidateADsToTest.drop(['FullProteinSeq'], axis=1)

    print("Saving output to: "+'../output/'+exportfilename)
    CandidateADsToTest.to_csv('../output/'+exportfilename)

    return CandidateADsToTest

# +
# Function to not combine the Predicted Tiles, just return tiles of predicted ADs--look for overlapping tiles
def TilesOfPredictedADs(inputfilename, Sequences_to_Test_DF, TEMPindx, exportfilename):
    print("Returning tiles that have not been aggregated!")
    CombinedPredictions = Sequences_to_Test_DF[TEMPindx]
    return CombinedPredictions
#     genenames, starts,ends, seqs,regionlengths,RegionType = [],[],[],[],[],[]
#     for gene in set(CombinedPredictions.GeneName):
#         indx = CombinedPredictions.GeneName == gene
#         tempDF = CombinedPredictions[indx]
#         tempDF = tempDF.reset_index()
#         display(tempDF)
#         start, end = min(tempDF.StartPosition)+1,max(tempDF.EndPosition)
#         starts.append(start)
#         ends.append(end)
#         regionlengths.append(end-start)
#         genenames.append(gene)
#         RegionType.append('Prediction')

#     # make a DF with geneNames and start and End positions of the overlapping windows
#     CandidateADsToTest = pd.DataFrame({'GeneName':genenames,'Start':starts,'End':ends,'Length':regionlengths,'RegionType':RegionType})

#     #add a column for full length AA seq
#     FullLengthProteinsTested = makeFullLengthProteinDF(inputfilename)

#     #print(FullLengthProteinsTested)
#     tempDict = dict(zip(FullLengthProteinsTested.GeneName,FullLengthProteinsTested.AAseq))
#     print(len(tempDict))
#     tempSeries = CandidateADsToTest.GeneName
#     print(tempSeries)
#     CandidateADsToTest['FullProteinSeq'] = tempSeries.map(tempDict)


#     ## pull out regions of these TFs
#     # CandidateADsToTest['Region2Test'] = CandidateADsToTest.apply(lambda row: row.FullProteinSeq[row.Start:row.End], axis=1)
#     #correct for counting starting at 0 in python and 1 in real world
#     ProteinRegionSeqs =[]
#     for i, Start, End, FullProteinSeq in zip(range(len(CandidateADsToTest.index)),CandidateADsToTest["Start"],CandidateADsToTest["End"],CandidateADsToTest["FullProteinSeq"]):
#         start =Start
#         end = End
#         Region = FullProteinSeq[start:end]
#         ProteinRegionSeqs.append(Region)
#     CandidateADsToTest['ProteinRegionSeq']=ProteinRegionSeqs
#     CandidateADsToTest = CandidateADsToTest.drop(['FullProteinSeq'], axis=1)

#     #print("Saving output to: "+'../output/'+exportfilename)
#     #CandidateADsToTest.to_csv('../output/'+exportfilename)

#     return CandidateADsToTest
# -

def make_predictions(folder_name="predictions/",
                        inputfilename="../data/LambertTFs.fasta", 
                        slope=1,
                        lower_corner_c="VP16",lower_corner_h="VP16",
                        upper_corner_c="CITED2",upper_corner_h="CITED2",
                        LowerCorner_slope1=0, LowerCorner_slope2=0,
                        UpperCorner_slope1=1, UpperCorner_slope2='inf',
                        composition=["W","F","Y","L"],
                        window_size=39,
                        window_spacing=1,
                        propset=["Charge","AllHydros"],
                        plot_new_red_points = True,
                        aggregate = True):
    
    exportfilename=return_exportfilename(folder_name,
                        inputfilename, 
                        slope,
                        lower_corner_c,lower_corner_h,
                        upper_corner_c,upper_corner_h,
                        LowerCorner_slope1, LowerCorner_slope2,
                        UpperCorner_slope1, UpperCorner_slope2,
                        composition,
                        window_size=39,
                        window_spacing=1,
                        propset=["Charge","AllHydros"])
                        
                        
    #exportfilename="testingmultipleslopes"

    lower_corner_c=get_bounds(lower_corner_c,composition)[0] * window_size / 39
    lower_corner_h=get_bounds(lower_corner_h,composition)[1] * window_size / 39
    
    upper_corner_c=get_bounds(upper_corner_c,composition)[0] * window_size / 39
    upper_corner_h=get_bounds(upper_corner_h,composition)[1] * window_size / 39
    
    Sequences_to_Test = makeTilingDF(inputfilename,window_size, window_spacing,composition)
    FullLengthProteinsTested = makeFullLengthProteinDF(inputfilename)
    PredictedTiles = maskproteome(inputfilename,composition, Sequences_to_Test, LowerCorner_slope1, LowerCorner_slope2,
                        UpperCorner_slope1, UpperCorner_slope2, slope, lineparameters=[[lower_corner_c,lower_corner_h],[upper_corner_c,upper_corner_h]],propset=["Charge","AllHydros"],
                        window_size = window_size, plot_new_red_points = plot_new_red_points)
    if aggregate == True:
        return AggregateTilesIntoPredictedADs(inputfilename, Sequences_to_Test,PredictedTiles,exportfilename)
    else:
        return TilesOfPredictedADs(inputfilename, Sequences_to_Test,PredictedTiles,exportfilename)

def get_bounds(AD_name,composition):

    if type(AD_name)==int:
        return [AD_name,AD_name]

    WFYL_GoldStandard=pd.read_csv("../data/Gold Standard AD List With Counts.csv")
    WFYL_GoldStandard["Charge"]=WFYL_GoldStandard.apply(lambda row: row['Sequence'].count('K')+row['Sequence'].count('R')-row['Sequence'].count('D')-row['Sequence'].count('E'),axis=1)
    # for aa in ["W","F","Y","L"]:
    #         WFYL_GoldStandard[aa]= WFYL_GoldStandard.apply(lambda row: row['Sequence'].count(aa),axis=1)
    # WFYL_GoldStandard['AllHydros'] = WFYL_GoldStandard[["W","F","Y","L"]].sum(axis=1)
    

    if(str(AD_name)=="VP16"):
        AD=WFYL_GoldStandard[(WFYL_GoldStandard["Charge"]==-13)&(WFYL_GoldStandard["AllHydros"]==7)&(WFYL_GoldStandard["Length"]==46)]
        AD=AD.drop(axis=1,labels=["Charge","AllHydros"])

    elif(str(AD_name)=="CITED2"):
        AD=WFYL_GoldStandard[(WFYL_GoldStandard["Charge"]==-9)&(WFYL_GoldStandard["AllHydros"]==10)&(WFYL_GoldStandard["GeneName"]==AD_name)]
        AD=AD.drop(axis=1,labels=["Charge","AllHydros"])
        
    else:
        print("ERROR: this function currently can't find bounds for AD other than VP16 or CITED2")
        return 0

    AD["Charge"]=AD.apply(lambda row: row['Sequence'].count('K')+row['Sequence'].count('R')-row['Sequence'].count('D')-row['Sequence'].count('E'),axis=1)
    AD_charge=AD["Charge"].iloc[0]

    AD["AllHydros"]=AD[composition].sum(axis=1)
    AD_AllHydros=AD["AllHydros"].iloc[0]
    
    return [AD_charge,AD_AllHydros]


def compare_to_random(Nrepeats=1,
    inputfilepath = "../data/LambertTFs.fasta",
    outputfilepath="../output/predictions/LambertTFs_s_1e-05_lcc_-40_lch_000_ucc_CITED2_uch_000_comp_WFYL_tl_039_ws_001",
    gs_list="../data/Gold Standard AD List.csv"):

    humanTFs = makeFullLengthProteinDF(inputfilepath)
    humanTFs['uniprotID']= humanTFs.apply(lambda row: row['GeneName'].split('|')[1],axis=1)

    if type(outputfilepath) == str:
        predictionDF=pd.read_csv(outputfilepath)
    else:
        predictionDF=outputfilepath
    if "uniprotID" not in predictionDF.columns:
        predictionDF['uniprotID']= predictionDF.apply(lambda row: row['GeneName'].split('|')[1],axis=1)

    CurrentUniprotADSet = pd.read_csv(gs_list)

    # find length distribution of predicted AD regions in the prediction set
    LengthSet = []
    for UPid in set(predictionDF.uniprotID):
        tempDF = predictionDF[predictionDF.uniprotID==UPid]
        templist =[]  
        for ProteinRegionSeq in tempDF["ProteinRegionSeq"]:
            templist.append(len(ProteinRegionSeq))
        LengthSet.append(templist)
    LengthDist = list(predictionDF['Length'].values)
    Nregions = len(LengthDist)
    NTFs2sample = len(set(predictionDF.uniprotID)) #when many TFs have multiple predictions it is more correct randomly sample multiple regions per TF 20200825
    HoldOverlapCounts = np.zeros((Nrepeats,Nregions))
    for j in np.arange(Nrepeats):
        SampledTFRegions =[]
        OverlapCounter = 0
        # randomly sample TFs
        TFsubset = humanTFs.sample(n=NTFs2sample)# May need to update this to be NTFs2sample 20200825
        TFsubset = TFsubset.reset_index()
        
        # cycle through TFs
        for i,uniprotID,AAseq in zip(np.arange(0,len(TFsubset.index)),TFsubset["uniprotID"],TFsubset["AAseq"]): # note I choose the random regions in the same order each time, but the TFs are randomized so it should be fine
            ## if TF has a known AD
            ## note, this only generates a random region of the TF if the selected TF has a known AD. Otherwise, the chance of overlap is Zero and to save time,I do not choose a random region
            if uniprotID in CurrentUniprotADSet.uniprotID.values:
                for length in LengthSet[i]: # for each length in the set of lengths for that TF
                # pick a random starting location--uniformly sample TF length
                # check that the region length is shorter than the TF length
                    if len(AAseq)>length:
                        randStart = random.sample(np.arange(len(AAseq)-length).tolist(),k=1)[0]
                        End = randStart+length    
                    else:
                        # if the expected region is longer than the sampled TF, use the full TF length
                        print('random region length is longer than this TF')
                        randStart, End = 0, len(AAseq)
                    # assess if sampled regions overlap know ADs
                    indx = CurrentUniprotADSet.uniprotID==uniprotID
                    TFwithKnownADs =  CurrentUniprotADSet[indx]
                    # Because a TF can have multiple AD entries, cycle through these entries
                    for KnownAD_Start,KnownAD_End in zip(TFwithKnownADs["Start"],TFwithKnownADs["End"]):
                        KnownStart = KnownAD_Start
                        KnownEnd = KnownAD_End
        #                 print randStart
        #                 print KnownStart
                        if (randStart<=KnownEnd)&(End>KnownStart):
                            HoldOverlapCounts[j,i]+=1
    # collapese hits from each TF to avoid double counting duplicates
    hits = np.sum(HoldOverlapCounts>0,axis=1)

    ## actual predictions
    predictionCounts,overlapcounter = np.zeros((1,Nregions)),0
    # overlap_GeneName,overlap_Sequence=[],[]
    for i,entry_uniprotID,entry_Start,entry_End,entry_GeneName in zip(np.arange(0,len(predictionDF.index)),predictionDF["uniprotID"],predictionDF["Start"],predictionDF["End"],predictionDF["GeneName"]): # note I choose the random regions in the same order each time, but the TFs are randomized so it should be fine    ## if TF has a known AD
        ## note, this only generates a random region of the TF if the selected TF has a known AD. Otherwise, the chance of overlap is Zero and to save time,I do not choose a random region
        if entry_uniprotID in CurrentUniprotADSet.uniprotID.values:
            predictedStart = entry_Start
            predictedEnd = entry_End
            # assess if sampled regions overlap know ADs
            indx = CurrentUniprotADSet.uniprotID==entry_uniprotID
            TFwithKnownADs =  CurrentUniprotADSet[indx]
            #print(TFwithKnownADs)
            # Because a TF can have multiple AD entries, cycle through these entries
            for KnownAD_Start,KnownAD_End,GeneName,Sequence in zip(TFwithKnownADs["Start"],TFwithKnownADs["End"],TFwithKnownADs["GeneName"],TFwithKnownADs["Sequence"]):
                KnownStart = KnownAD_Start
                KnownEnd = KnownAD_End
                print(entry_GeneName)
                print('%i - %i Known'%( KnownStart,KnownEnd))
                print('%i - %i Predicted'%( predictedStart,predictedEnd)) 
                if (predictedStart<=KnownEnd)&(predictedEnd>KnownStart):
                    overlapcounter +=1
                    predictionCounts[0,i]+=1
                    #overlap_GeneName.append(GeneName)
                    #overlap_Sequence.append(Sequence)
    Npredictions =np.sum(predictionCounts>0,axis=1)
    print('There are %i predictions total. %i overlap known ADs'%(Nregions,Npredictions))
    # print(overlap_GeneName)
    # print(overlap_Sequence)

    figure(dpi=200)

    bins = np.arange(-.5,Npredictions+5.5,1)
    ys =plt.hist(hits,bins,density=1)
    plt.plot([Npredictions,Npredictions],[0,.4])


    plt.xlabel('Number of overlaps')
    plt.ylabel('Density')
    if type(outputfilepath) == str:
        plt.savefig('../output/compare_to_random/'+outputfilepath.split("/")[-1]+"_Nrepeats_"+str(Nrepeats)+"_comp_to_"+gs_list.split("/")[-1].split(".")[0]+".png")
    plt.show()
    sum(hits>Npredictions)
    print('There are %i predictions total. %i overlap known ADs'%(Nregions,Npredictions)) 
    print('Ran %i permutations'%Nrepeats) 
    print('Greatest number of times random regions overlap known ADs is %s'% max(hits)) 
    print('nominal p value is %f'% (1.0*sum(hits>Npredictions)/ Nrepeats)) 

    return pd.DataFrame({"NumberOfPredictedADs":[Nregions],
                        "NumberPredicted_GSOverlaps":Npredictions,
                        "ProportionPredictedAD_GoldStandardOverlaps":Npredictions/Nregions,
                        "MaxNumberRandom_GSOverlaps":max(hits),
                        "p_value":[1.0*sum(hits>Npredictions)/ Nrepeats]})

# make_predictions_and_compare_to_random():
#     make_predictions(inputfilename="'../data/LambertTFs.fasta", slope=1,lower_corner_c=-13,lower_corner_h=7,upper_corner_c=-9,upper_corner_h=10,composition=["W","F","Y","L"],window_size=39,window_spacing=1)


