# VIRANET TOOLS #
###############################################################################
# objSto
# Strips stoichiometric coefficients from the objectives reactions
# usage
# Inputs:
# HVM               Integrated host-virus model
# HostRxn           Host objective reaction, either:
#                   - Index value of reaction in Model.reactions        [int]
#                   - Reaction ID of the host-objective reaction        [str]
#
# Outputs:
# hostSto           Stoichiometric coefficients with metabolites for host objective
# virusSto          Stoichiometric coefficients with metabolites for virus objective
#

def objSto(HVM,HostRxn):
    "Differential usage of amino acids and nucleotides"
    # [1] Initial Setup
    # Function Dependencies
    import pandas   as pd
    import numpy    as np
    # Identify the host objective reaction
    try:
        intTest = int(HostRxn)
        hostIdx = HostRxn
    except:
        for ii in range(len(HVM.reactions)):
            if HostRxn in str(HVM.reactions[ii]):
                hostIdx = ii
    virusIdx    = len(HVM.reactions) - 1    # ViraNet(c) appends virus reaction to end in genHVM.py
    objIdx      = [hostIdx,virusIdx]
    hostID      = HVM.reactions[hostIdx].id
    virusID     = HVM.reactions[virusIdx].id
    # Condition: ensure virus reaction is virus objective
    if '_prodrxn_VN' not in virusID:
        raise ValueError('Unsupported objective, unable to analyse: refer to README')
    # Convert model into an array
    m       = HVM.to_array_based_model()
    # Create data frame
    mFrame  = pd.DataFrame(
        data    = m.S.todense(),
        columns = m.reactions.list_attr("id"),
        index   = m.metabolites.list_attr("id")
    )
    # [2] Stoichiometric Data
    # Strip the host and virus objective function stoichiometric coefficients
    hostS   = mFrame[hostID]
    virusS  = mFrame[virusID]
    # Record only non zero values
    hostSto     = hostS.loc[~(hostS==0)]
    virusSto    = virusS.loc[~(virusS==0)]
    # [3] Create output
    hostSto     = pd.DataFrame(hostSto)
    virusSto    = pd.DataFrame(virusSto)
    return (hostSto,virusSto)
###############################################################################

###############################################################################
# boundCheck
# Checks a model for any arbitrarily large bounds, removes these 
# and replaces with infinite bounds
# Work based upon Kelk et al 2012
# Kelk, S. M., Olivier, B. G., Stougie, L., & Bruggeman, F. J. (2012). Optimal flux spaces of genome-scale stoichiometric models are determined by a few subnetworks. Scientific Reports, 2, 580. http://doi.org/10.1038/srep00580
# Inputs:
# Model         User-supplied model file   [.mat,.xml]
#
# Outputs:
# infModel      Model with bound corrections applied
def boundCheck(Model):
    "Reaction bound checker"
    # [1] Initial Setup
    # Function Dependencies
    import numpy    as np
    # Create pointer
    altModel = Model
    # [2] Identify and correct arbitrarily large reaction bounds
    for ii in range(len(altModel.reactions)):
        # Temporarily record the lower and upper bounds
        tmpLb   = altModel.reactions[ii].lower_bound
        tmpUb   = altModel.reactions[ii].upper_bound
        # Conditional statement
        if tmpLb <= -1000:
            altLb   = -np.inf
            altModel.reactions[ii].lower_bound = altLb
        if tmpUb >= 1000:
            altUb   = np.inf
            altModel.reactions[ii].upper_bound = altUb
    # [3] Output model
    return altModel
###############################################################################
# rangeCalculator
# Calculates , using FVA results for host and virus, the flux range to use in the 
# host-derived enforcement analysis
# Inputs:
# hostIdx           Index (model.reactions) for the host-objective reaction
# virusIdx          Index (model.reactions) for the virus-objective reaction
# Optional Inputs
# solver            Declare solver to use for cobrapy: default is cglpk
# Outputs
# enfVirus          Vector of virus optima values with additional host-constraint
# maxEnfBound       Maximum enf bound
# minEnfBound       Minimum enf bound
def rangeCalculator(HVM,hostIdx,virusIdx,solver):
    "Enforcement bound creator"
    # [1] Initial Setup
    # Function Dependencies
    import cobra
    import numpy    as np
    import pandas   as pd
    # Create the virus optima vector
    enfVirus    = np.zeros((len(HVM.reactions),1))
    maxEnfBound = np.zeros((len(HVM.reactions),1))
    minEnfBound = np.zeros((len(HVM.reactions),1))
    # [2] Perform FVA for each objective
    # Objective reactions
    hostObj     = HVM.reactions[hostIdx]
    virusObj    = HVM.reactions[virusIdx]
    # Host Optimisation
    HVM.change_objective(hostObj)
    # Ensure no flux can go through the virus reaction
    # Store the bounds
    virusLb     = HVM.reactions[virusIdx].lower_bound
    virusUb     = HVM.reactions[virusIdx].upper_bound
    HVM.reactions[virusIdx].lower_bound = 0
    HVM.reactions[virusIdx].upper_bound = 0
    # FVA
    varHost     = cobra.flux_analysis.flux_variability_analysis(HVM,solver=solver)
    # Return virus objective bounds
    HVM.reactions[virusIdx].lower_bound = virusLb
    HVM.reactions[virusIdx].upper_bound = virusUb
    # Virus Optimisation
    HVM.change_objective(virusObj)
    # Ensure no flux can go through the host reaction
    # Store the bounds
    hostLb      = HVM.reactions[hostIdx].lower_bound
    hostUb      = HVM.reactions[hostIdx].upper_bound
    HVM.reactions[hostIdx].lower_bound = 0
    HVM.reactions[hostIdx].upper_bound = 0
    # FVA
    varVirus    = cobra.flux_analysis.flux_variability_analysis(HVM,solver=solver)
    # Return host objective bounds
    HVM.reactions[hostIdx].lower_bound = hostLb
    HVM.reactions[hostIdx].upper_bound = hostUb
    # Create data frames
    hostFVA     = pd.DataFrame.from_dict(varHost)
    hostFVA     = hostFVA.transpose()
    virusFVA    = pd.DataFrame.from_dict(varVirus)
    virusFVA    = virusFVA.transpose()
    # [3] Condition statements to determine the calculation and FVA steps    
    # Initiate loop
    for ii in range(len(HVM.reactions)):
        # Create temporary host and virus max|min variables
        hostMax     = varHost[HVM.reactions[ii].id]['maximum']
        hostMin     = varHost[HVM.reactions[ii].id]['minimum']
        virusMax    = varVirus[HVM.reactions[ii].id]['maximum']
        virusMin    = varVirus[HVM.reactions[ii].id]['minimum']
        # Record the upper and lower bounds
        tmpLb       = HVM.reactions[ii].lower_bound
        tmpUb       = HVM.reactions[ii].upper_bound
        # Conditional: H+ > V+ && H- < V-
        if (hostMax > virusMax) and (hostMin < virusMin):
            #####################################################################
            # Calculation of bounds for condition [1]
            enfMax1 = hostMax
            enfMin1 = (hostMax - ((hostMax - virusMax) / 2))
            # Apply bounds to reaction
            HVM.reactions[ii].lower_bound   = enfMin1
            HVM.reactions[ii].upper_bound   = enfMax1
            # Zero-bound host
            HVM.reactions[hostIdx].lower_bound = 0
            HVM.reactions[hostIdx].upper_bound = 0
            # Optimize for virus
            HVM.change_objective(virusObj)
            sol = HVM.optimize(objective_sense='maximize',solver=solver)
            # Record the optima
            zMax    = sol.f
            # Return host bounds to original
            HVM.reactions[hostIdx].lower_bound = hostLb
            HVM.reactions[hostIdx].upper_bound = hostUb
            # Return reaction bounds to original
            HVM.reactions[ii].lower_bound   = tmpLb
            HVM.reactions[ii].upper_bound   = tmpUb
            # Calculation of bounds for condition [2]
            enfMin2 = hostMin
            enfMax2 = (hostMin - ((hostMin - virusMin) / 2))
            # Apply bounds to reaction
            HVM.reactions[ii].lower_bound   = enfMin2
            HVM.reactions[ii].upper_bound   = enfMax2
            # Zero-bound host
            HVM.reactions[hostIdx].lower_bound = 0
            HVM.reactions[hostIdx].upper_bound = 0
            # Optimize for virus
            HVM.change_objective(virusObj)
            sol = HVM.optimize(objective_sense='maximize',solver=solver)
            # Record the optima
            zMin    = sol.f
            # Return host bounds to original
            HVM.reactions[hostIdx].lower_bound = hostLb
            HVM.reactions[hostIdx].upper_bound = hostUb
            # Return reaction bounds to original
            HVM.reactions[ii].lower_bound   = tmpLb
            HVM.reactions[ii].upper_bound   = tmpUb
            # COMPARISON #
            # Compare zMax and zMin to find which is smallest
            if zMax < zMin:
                enfVirus[ii] = zMax
                # Record the bound
                maxEnfBound[ii] = enfMax1
                minEnfBound[ii] = enfMin1
            elif zMin < zMax:
                enfVirus[ii] = zMin
                # Record the bound
                maxEnfBound[ii] = enfMax2
                minEnfBound[ii] = enfMin2
            else:
                enfVirus[ii] = zMax
                # Record the bound
                maxEnfBound[ii] = enfMax1
                minEnfBound[ii] = enfMin1
        # Conditional: H+ > V+
        else:
            if hostMax > virusMax:
                #####################################################################
                # Calculation of bounds
                enfMax  = hostMax
                enfMin  = (hostMax - ((hostMax - virusMax) / 2))
                # Apply bounds to reaction
                HVM.reactions[ii].lower_bound   = enfMin
                HVM.reactions[ii].upper_bound   = enfMax
                # Zero-bound host
                HVM.reactions[hostIdx].lower_bound = 0
                HVM.reactions[hostIdx].upper_bound = 0
                # Optimize for virus
                HVM.change_objective(virusObj)
                sol = HVM.optimize(objective_sense='maximize',solver=solver)
                # Record the optima
                enfVirus[ii]    = sol.f
                # Return host bounds to original
                HVM.reactions[hostIdx].lower_bound = hostLb
                HVM.reactions[hostIdx].upper_bound = hostUb
                # Return reaction bounds to original
                HVM.reactions[ii].lower_bound   = tmpLb
                HVM.reactions[ii].upper_bound   = tmpUb
                # Record the bound
                maxEnfBound[ii] = enfMax
                minEnfBound[ii] = enfMin
                #####################################################################
            else:
                # Conditional: H- < V-
                if hostMin < virusMin:
                    #####################################################################
                    enfMin  = hostMin
                    enfMax  = (hostMin - ((hostMin - virusMin) / 2))
                    # Apply bounds to reaction
                    HVM.reactions[ii].lower_bound   = enfMin
                    HVM.reactions[ii].upper_bound   = enfMax
                    # Zero-bound host
                    HVM.reactions[hostIdx].lower_bound = 0
                    HVM.reactions[hostIdx].upper_bound = 0
                    # Optimize for virus
                    HVM.change_objective(virusObj)
                    sol = HVM.optimize(objective_sense='maximize',solver=solver)
                    # Record the optima
                    enfVirus[ii]    = sol.f
                    # Return host bounds to original
                    HVM.reactions[hostIdx].lower_bound = hostLb
                    HVM.reactions[hostIdx].upper_bound = hostUb
                    # Return reaction bounds to original
                    HVM.reactions[ii].lower_bound   = tmpLb
                    HVM.reactions[ii].upper_bound   = tmpUb
                    # Record the bound
                    maxEnfBound[ii] = enfMax
                    minEnfBound[ii] = enfMin
                    #####################################################################
                else:
                    #####################################################################
                    enfMax  = hostMax
                    enfMin  = hostMin
                    # Apply bounds to reaction
                    HVM.reactions[ii].lower_bound   = enfMin
                    HVM.reactions[ii].upper_bound   = enfMax
                    # Zero-bound host
                    HVM.reactions[hostIdx].lower_bound = 0
                    HVM.reactions[hostIdx].upper_bound = 0
                    # Optimize for virus
                    HVM.change_objective(virusObj)
                    sol = HVM.optimize(objective_sense='maximize',solver=solver)
                    # Record the optima
                    enfVirus[ii]    = sol.f
                    # Return host bounds to original
                    HVM.reactions[hostIdx].lower_bound = hostLb
                    HVM.reactions[hostIdx].upper_bound = hostUb
                    # Return reaction bounds to original
                    HVM.reactions[ii].lower_bound   = tmpLb
                    HVM.reactions[ii].upper_bound   = tmpUb
                    # Record the bound
                    maxEnfBound[ii] = enfMax
                    minEnfBound[ii] = enfMin
                    #####################################################################
    # [4] Output
    return (enfVirus,maxEnfBound,minEnfBound)
###############################################################################