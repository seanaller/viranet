###############################################################################
# main.py performs analysis on a user-defined host model, and a user-
# supplied virus genome file (from NCBI > Genbank)
# A host-virus integrated model is created and the following analysis performed:
# - Comparison of flux distribution
# - Analysis of the effect of single-reaction knockouts on virus optima
# - Analysis of the effect of single-reaction host-derived flux enforcement on
#   virus optima
###############################################################################
# Sean Aller: 2016 05 23
# S-D.Aller@warwick.ac.uk
###############################################################################
# Setup of workspace
# Python Dependencies
import  cobra
import  numpy       as np
import  pandas      as pd
import  os
import  re
import  csv
import  time
# Intra-package dependencies
from    viranet     import generation
from    viranet     import analysis
from    viranet     import info
from    viranet     import tools

###############################################################################
# singleAnalysis
# Performs analysis on a integrated host-virus model, for a user-defined
# metabolic model and user-defined virus
# Inputs:
# Model             User-supplied model file for desired host   [.mat,.xml]
# VirusGB           User-supplied GenBank file (NCBI) for desired virus [.gb, .txt]
# HostRxn           Host objective reaction, either:
#                   - Index value of reaction in Model.reactions        [int]
#                   - Reaction ID of the host-objective reaction        [str]
# Optional Inputs:
# solver            User-specified solver (default = CGLPK)
# preanalysis       User-specified pre-analysis of model (default = False)
#                   - Checking for arbitrarily large (> 1000) reaction flux bounds
# usefva            Use the Host FBA result (False) or use the FVA range as
#                   the additional constraint (True)
#
# Outputs:
# ViraNet_Analysis  Full analysis of the HVM, written to .csv files

def singleAnalysis(Model,VirusGB,HostRxn,solver="cglpk",preanalysis=False,usefva=False,userange=False):
    "Analysis of Host-Virus model"
    # Generate the neccessary virus objective function and integrated host-virus model
    # require for the analysis
    # [1] VBOF Generation
    virusReaction = generation.VBOF(VirusGB)
    print("[1] VBOF Generated")

    # [2] Intergrated Host-Virus Model (HVM)
    (virusModel) = generation.HVM(Model,virusReaction)
    print("[2] HVM Generated")

    ###################################################
    ## PUBLICATION SPECIFIC REMOVE BEFORE SUBMISSION ##
    CTPS2   = virusModel.reactions[1084]
    virusModel.remove_reactions(CTPS2)
    ## PUBLICATION SPECIFIC REMOVE BEFORE SUBMISSION ##
    ###################################################

    # OPTIONAL STEPS #
    if preanalysis == True:
        # [I] Checking model for arbitrarily large bounds
        virusModel = tools.boundCheck(virusModel)
        print("[2A] HVM Reaction Bounds Checked")
    
    # [3] Host-Virus Optimisation
    # Optimise the original model
    (objIdx,hostF,hostX,virusF,virusX) = analysis.optimise(virusModel,HostRxn,solver)
    print("[3] HVM Optimisation Complete")

    # [4] Host-Virus Differential Analysis
    (rel_dif, abs_dif) = analysis.differential(virusModel,HostRxn)
    print("[4] HVM Differential Analysis Complete")

    # [5] Host-Virus Comparison
    (hvmComp,hvmStat,hostXd,virusXd) = analysis.compare(objIdx,hostX,virusX)
    print("[5] HVM Comparison Complete")

    # [6] Host-Virus Variability Analysis
    (hostFVA,virusFVA) = analysis.variability(virusModel,HostRxn,solver)
    print("[6] HVM Variability Analysis Complete")

    # [7] Knockout Analysis
    (koVirus)   = analysis.knockout(virusModel,HostRxn,solver)
    print("[7] HVM Knockout Analysis Complete")

    # [8] Enforcement Analysis
    # Conditional on userange
    if userange == False:
        (enfVirus, _, _) = analysis.enforce(virusModel,hostX,HostRxn,solver,usefva,userange)
    elif userange == True:
        (enfVirus,maxEnfBound,minEnfBound) = analysis.enforce(virusModel,hostX,HostRxn,solver,usefva,userange)
    print("[8] HVM Host-Derived Enforcement Analysis Complete")

    # [9] Function output: setup
    # Identify the host objective reaction
    try:
        intTest = int(HostRxn)
        hostIdx = HostRxn
    except:
        for ii in range(len(virusModel.reactions)):
            if HostRxn in str(virusModel.reactions[ii]):
                hostIdx = ii
    # Record the other reactions
    virusIdx    = len(virusModel.reactions) - 1    # ViraNet(c) appends virus reaction to end in genHVM.py
    objIdx      = [hostIdx,virusIdx]
    hostID      = virusModel.reactions[hostIdx].id
    virusID     = virusModel.reactions[virusIdx].id
    # Condition: ensure virus reaction is virus objective
    if '_prodrxn_VN' not in virusID:
        raise ValueError('Unsupported objective, unable to analyse: refer to README')
    # Create array-model and obtain neccessary information
    m       = virusModel.to_array_based_model()
    mRxns   = m.reactions.list_attr("id")
    # Remove the objective reactions
    if hostID in mRxns:
        mRxns.remove(hostID)
    if virusID in mRxns:
        mRxns.remove(virusID)
    # Obtain subsystem information
    subSystems = []
    for ii in range(len(virusModel.reactions)):
        if not(ii in objIdx):
            subSystems.append(virusModel.reactions[ii].subsystem)

    # [10] singleAnalysis Output Generation
    # Virus name and date
    date        = (time.strftime("%d_%m_%Y"))
    idVirus     = virusID.replace("_prodrxn_VN","")
    idVirus     = "_" + idVirus + "_" + date

    # OUTPUT [1]: Host-Virus Optimisation Comparisons
    # Store as list
    hostFlux        = list(hvmComp[:,0])
    virusFlux       = list(hvmComp[:,1])
    hostSolX        = list(hostXd)
    virusSolX       = list(virusXd)
    # Store in dataframe
    hvmComparison   = pd.DataFrame([mRxns,subSystems,hostSolX,virusSolX,hostFlux,virusFlux])
    hvmComparison   = hvmComparison.transpose()
    hvmComparison.columns = ['Reaction','Subsystem','Host_Flux_mmol/gDW/h-1','Virus_Flux__mmol/gDW/h-1','Host_Flux_%','Virus_Flux_%']
    title = "ViraNet_Comparison" + idVirus + ".csv"
    # Write to csv
    hvmComparison.to_csv(title,index=False)

    # OUTPUT [2]: Host-Virus Optimisation Comparison Stats
    # Store as list
    compStats       = list(hvmStat)
    nameStats       = ['Upregulated','Downregulated','Activated','Inactivated','Reversed']
    # Store in dataframe
    hvmComparisonStats  = pd.DataFrame([nameStats,compStats])
    hvmComparisonStats  = hvmComparisonStats.transpose()
    hvmComparisonStats.columns  = ['Regulation_state','Num_Rxns']
    title = "ViraNet_ComparisonStats" + idVirus + ".csv"
    # Write to csv
    hvmComparisonStats.to_csv(title,index=False)

    # OUTPUT[3]: Host-virus variability analysis using FVA
    # Create title host
    title = "ViraNet_HostFVA" + idVirus + ".csv"
    # Write to csv
    hostFVA.to_csv(title,index=True)
    # Create title virus
    title = "ViraNet_VirusFVA" + idVirus + ".csv"
    # Write to csv
    virusFVA.to_csv(title,index=True)

    # OUTPUT[4]: Host-virus differential usage of amino acids and nucleotides
    # Create lists for the relative and absolute measurements
    relativeDiff    = list(rel_dif.values())
    absoluteDiff    = list(abs_dif.values())
    orderedMets     = list(rel_dif.keys())
    # Store in dataframe
    hvmDifferential     = pd.DataFrame([orderedMets,relativeDiff,absoluteDiff])
    hvmDifferential     = hvmDifferential.transpose()
    hvmDifferential.columns     = ['Obj_Met','Relative_d','Absolute_d']
    title = "ViraNet_Differential_Usage" + idVirus + ".csv"
    # Write to csv
    hvmDifferential.to_csv(title,index=False)

    # OUTPUT [5]: Virus optima and effects of knockout and host-derived enforcement
    # Store as list whilst removing objective reactions
    resultKO        = list(np.delete(koVirus,objIdx))
    resultEF        = list(np.delete(enfVirus,objIdx))
    if userange == True:
        resultMaxEnfBound   = list(np.delete(maxEnfBound,objIdx))
        resultMinEnfBound   = list(np.delete(minEnfBound,objIdx))
    # Convert to percentage of wild-type [unconstrained] virus optima
    for ii in range(len(resultKO)):
        if np.isnan(resultKO[ii]) == False:
            resultKO[ii]    = (resultKO[ii] / virusF) * 100
        elif np.isnan(resultKO[ii]) == True:
            resultKO[ii]    = 0
        elif resultKO[ii] is None:
            resultKO[ii]    = 0
        if np.isnan(resultEF[ii]) == False:
            resultEF[ii]    = (resultEF[ii] / virusF) * 100
        elif np.isnan(resultEF[ii]) == True:
            resultEF[ii]    = 0
        elif resultEF[ii] is None:
            resultEF[ii]    = 0
    # Store in dataframe
    if userange == False:
        hvmKoEnf        = pd.DataFrame([mRxns,subSystems,resultKO,resultEF])
        hvmKoEnf        = hvmKoEnf.transpose()
        hvmKoEnf.columns    = ['Reaction','Subsystem','Knockout_Optima(%_WT)','Enforcement_Optima(%_WT)']
        title = "ViraNet_Knockout_Enforcement" + idVirus + ".csv"
        # Write to csv
        hvmKoEnf.to_csv(title,index=False)
        #return (hvmComp, hvmStat, koVirus, enfVirus)
    elif userange == True:
        hvmKoEnf        = pd.DataFrame([mRxns,subSystems,resultKO,resultEF,resultMaxEnfBound,resultMinEnfBound])
        hvmKoEnf        = hvmKoEnf.transpose()
        hvmKoEnf.columns    = ['Reaction','Subsystem','Knockout_Optima(%_WT)','Enforcement_Optima(%_WT)','MaximumBound_Enf','MinimumBound_Enf']
        title = "ViraNet_Knockout_Enforcement" + idVirus + ".csv"
        # Write to csv
        hvmKoEnf.to_csv(title,index=False)
        #return (hvmComp, hvmStat, koVirus, enfVirus)
    print("ViraNet Analysis Complete")
###############################################################################

###############################################################################
# def multi(Model,FileArray,HostRxn):
#     "Analysis of multiple host-virus pairs"
###############################################################################
