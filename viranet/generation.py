###############################################################################
# generation.py generates the virus biomass objective function (VBOF), and the
# host-virus integrated model (HVM) for a given virus genome file (virusGB) and
# a given host metabolic model (hostModel) file (.mat,.xml)
###############################################################################
# Sean Aller: 2016 05 23
# S-D.Aller@warwick.ac.uk
###############################################################################
# Setup of workspace
import  cobra
from    cobra                   import Model, Reaction, Metabolite
import  numpy                   as np
from    viranet.info            import metDict, ntpsDict, aaDict, miscDict, ntpsMets, aaMets, N_A, k_atp, k_ppi
import  re

# Definitions
# Virus locations
# This indicate the 6 starting amino acids of the non-structural polyprotein
# for the supported flavivurses. These are used to determine the seperation
# of structural and non-structural polyproteins
DENVloc = 'DSGCVV'
ZIKVloc = 'DVGCSV'
# Metabolite definitions for final objective reaction creation
# This creates a mapping between the python-cobrapy and the ViraNet function
atp_c = Metabolite(metDict['atp'])
ctp_c = Metabolite(metDict['ctp'])
gtp_c = Metabolite(metDict['gtp'])
utp_c = Metabolite(metDict['utp'])
ala_c = Metabolite(metDict['A'])
arg_c = Metabolite(metDict['R'])
asn_c = Metabolite(metDict['N'])
asp_c = Metabolite(metDict['D'])
cys_c = Metabolite(metDict['C'])
gln_c = Metabolite(metDict['Q'])
glu_c = Metabolite(metDict['E'])
gly_c = Metabolite(metDict['G'])
his_c = Metabolite(metDict['H'])
ile_c = Metabolite(metDict['I'])
leu_c = Metabolite(metDict['L'])
lys_c = Metabolite(metDict['K'])
met_c = Metabolite(metDict['M'])
phe_c = Metabolite(metDict['F'])
pro_c = Metabolite(metDict['P'])
ser_c = Metabolite(metDict['S'])
thr_c = Metabolite(metDict['T'])
trp_c = Metabolite(metDict['W'])
tyr_c = Metabolite(metDict['Y'])
val_c = Metabolite(metDict['V'])
adp_c = Metabolite(metDict['adp'])
h2o_c = Metabolite(metDict['h2o'])
h_c   = Metabolite(metDict['h'])
pi_c  = Metabolite(metDict['Pi'])
ppi_c = Metabolite(metDict['PPi'])
###############################################################################
# Generation of the virus biomass objective function
# Inputs:
# VirusGB   User-supplied GenBank file (NCBI) for desired virus [.gb, .txt]

# Outputs:
# VBOF      Virus biomass objective function for desired virus

def VBOF(VirusGB):
    "Generate_VBOF"

    # [1] Initial Setup
    # Open virus file and parse contents
    with open(VirusGB,'rU') as vf:
        virusFile     = vf.readlines()
    for ii in range(len(virusFile)):
        virusFile[ii] = virusFile[ii].rstrip('\n')
    # Identify virus genera and define virusMethod
    # virusMethod denotes the method of VBOF generation to use, based upon the
    # virus genera
    virusData = str(virusFile)
    if "alphavirus".lower()     in virusData.lower():
        virusMethod = 1
    elif "flavivirus".lower()   in virusData.lower():
        virusMethod = 2
    else:
        raise ValueError('Unsupported virus, unable to create VBOF: refer to README')
    # VBOF Parameters: virusMethod dependent
    # Alphavirus
    if virusMethod == 1:
        # Copy number for viral genome                      [Cg]                # Source: Strauss, J. H., & Strauss, E. G. (1994).
        Cg  = 1
        # Copy number for viral structural polyprotein      [Csp]               # Source: Strauss, J. H., & Strauss, E. G. (1994).
        Csp = 240
        # Copy number for viral nonstructural polyprotein   [Cnp]               # Source: Strauss, J. H., & Strauss, E. G. (1994).
        Cnp = 1
    # Flavivirus
    elif virusMethod == 2:
        # Copy number for viral genome                      [Cg]                # Source: Mukhopadhyay, S., Kuhn, R. J., & Rossmann, M. G. (2005).
        Cg  = 1
        # Copy number for viral structural polyprotein      [Csp]               # Source: Mukhopadhyay, S., Kuhn, R. J., & Rossmann, M. G. (2005).
        Csp = 180
        # Copy number for viral nonstructural polyprotein   [Cnp]               # Source: Mukhopadhyay, S., Kuhn, R. J., & Rossmann, M. G. (2005).
        Cnp = 1
    # Virus Name
    # FUTURE UPDATE: Definition via blatimore classification and genera / species
    if "chikungunya".lower() in virusData.lower():
        virusName = "CHIKV"
        virusFull = "Chikungunya Virus"
    elif "semliki".lower() in virusData.lower():
        virusName = "SFV"
        virusFull = "Semliki Forest Virus"
    elif "sindbis".lower() in virusData.lower():
        virusName = "SINV"
        virusFull = "Sindbis Virus"
    elif "dengue".lower() in virusData.lower():
        virusName = "DENV"
        virusFull = "Dengue Virus"
        flavMeth  = 1                                                           # Indicator for nonstructural start
    elif "zika".lower() in virusData.lower():
        virusName = "ZIKV"
        virusFull = "Zika Virus"
        flavMeth  = 2                                                           # Which Virus Location paramter to use
    elif "eastern".lower() in virusData.lower():
        virusName = "EEEV"
        virusFull = "Eastern Equine Encephalitis Virus"
    elif "western".lower() in virusData.lower():
        virusName = "WEEV"
        virusFull = "Western Equine Encephalitis Virus"
    elif "venezuelan".lower() in virusData.lower():
        virusName = "VEEV"
        virusFull = "Venezuelan Equine Encephalitis Virus"
    else:
        raise ValueError('Unsupported virus, unable to create VBOF: refer to README')

    # [2] Sequence Identification
    # Genome Sequence: initial step is to identify start/end positions in the virus file
    startG  = [jj for jj, s in enumerate(virusFile) if 'ORIGIN' in s]
    endG    = [jj for jj, s in enumerate(virusFile) if '//' in s]
    startG  = int(''.join(map(str,startG)))
    endG    = int(''.join(map(str,endG)))
    startG  = startG + 1
    # Store genome sequence
    regex           = re.compile('[^a-zA-Z]')
    virusGenome     = str(''.join(virusFile[startG:endG]))
    virusGenome     = regex.sub('',virusGenome)
    # Polyprotein sequences: split into structural and nonstructural
    # This step is virus genera dependent (supported viruses only, see README)
    # Alphavirus
    if virusMethod == 1:
        # Structural polyprotein identification
        tempS1          = [jj for jj, s in enumerate(virusFile) if '/product="structural polyprotein"' in s]
        tempStruct      = virusFile[tempS1[0]:]
        tempS2          = [jj for jj, s in enumerate(tempStruct) if '/translation' in s]
        tempStruct      = tempStruct[tempS2[0]:]
        tempS3          = [jj for jj, s in enumerate(tempStruct) if 'gene' in s]
        tempStruct      = tempStruct[:tempS3[0]]
        # Clean-up and Store
        structReg       = re.compile('/translation=')
        virusStruct     = str(tempStruct)
        virusStruct     = structReg.sub('',virusStruct)
        virusStruct     = regex.sub('',virusStruct)
        # Nonstructural polyprotein identification
        tempN1          = [jj for jj, s in enumerate(virusFile) if '/product="nonstructural polyprotein"' in s]
        tempNP          = virusFile[tempN1[0]:]
        tempN2          = [jj for jj, s in enumerate(tempNP) if '/translation' in s]
        tempNP          = tempNP[tempN2[0]:]
        tempN3          = [jj for jj, s in enumerate(tempNP) if 'gene' in s]
        tempNP          = tempNP[:tempN3[0]]
        # Clean-up and Store
        npReg           = re.compile('/translation=')
        virusNonStruct  = str(tempNP)
        virusNonStruct  = npReg.sub('',virusNonStruct)
        virusNonStruct  = regex.sub('',virusNonStruct)
    # Flavivirus
    elif virusMethod == 2:
        # Polyprotein identification
        tempP1          = [jj for jj, s in enumerate(virusFile) if '/product="flavivirus polyprotein"' in s]
        tempPol         = virusFile[tempP1[0]:]
        tempP2          = [jj for jj, s in enumerate(tempPol) if '/translation' in s]
        tempPol         = tempPol[tempP2[0]:]
        tempP3          = [jj for jj, s in enumerate(tempPol) if 'gene' in s]
        tempPol         = tempPol[:tempP3[0]-1]
        npReg           = re.compile('/translation=')
        virusPoly       = str(tempPol)
        virusPoly       = npReg.sub('',virusPoly)
        virusPoly       = regex.sub('',virusPoly)
        # Structural and nonstructural polyprotein identification
        # flavMeth conditional variable: [1] Dengue virus; [2] Zika virus
        if flavMeth == 1:
            nsInd           = re.search(DENVloc, virusPoly).start()
            virusNonStruct  = virusPoly[nsInd:]
            virusStruct     = virusPoly[:nsInd]
        elif flavMeth == 2:
            nsInd           = re.search(ZIKVloc, virusPoly).start()
            virusNonStruct  = virusPoly[nsInd:]
            virusStruct     = virusPoly[:nsInd]
    # No supported virus detected
    else:
        raise ValueError('Unsupported virus, unable to create VBOF: refer to README')

    # [3] Precursor frequency
    # Genome                            [Nucleotides]
    countA  = virusGenome.count('a')
    countC  = virusGenome.count('c')
    countG  = virusGenome.count('g')
    countU  = virusGenome.count('t')    # Base 'T' is psuedo for base 'U'
    antiA   = countU
    antiC   = countG
    antiG   = countC
    antiU   = countA
    # Structural polyprotein            [Amino Acids]
    structCount     = np.zeros((20,1))
    for ii in range(len(aaMets)):
        structCount[ii,0] = virusStruct.count(aaMets[ii])
    # Nonstructural polyprotein         [Amino Acids]
    nonstructCount  = np.zeros((20,1))
    for ii in range(len(aaMets)):
        nonstructCount[ii,0] = virusNonStruct.count(aaMets[ii])
    # Count summation
    totNTPS     = (Cg * (countA + countC + countG + countU + antiA + antiC + antiG + antiU))
    totAA       = (structCount * Csp) + (nonstructCount * Cnp)

    # [4] VBOF Calculations
    # Nucleotides
    # mol.ntps/mol.virus
    V_a = (Cg*(countA + antiA))
    V_c = (Cg*(countC + antiC))
    V_g = (Cg*(countG + antiG))
    V_u = (Cg*(countU + antiU))
    # g.ntps/mol.virus
    G_a = V_a * ntpsDict["atp"]
    G_c = V_c * ntpsDict["ctp"]
    G_g = V_g * ntpsDict["gtp"]
    G_u = V_u * ntpsDict["utp"]
    # Amino Acids
    # mol.aa/mol.virus
    V_aa    = np.zeros((20,1))
    for ii in range(len(aaMets)):
        V_aa[ii,0] = totAA[ii]
    # g.a/mol.virus
    G_aa    = np.zeros((20,1))
    for ii in range(len(aaMets)):
        G_aa[ii,0] = V_aa[ii] * aaDict[aaMets[ii]]
    # Total genomic and proteomic molar mass
    M_v     = (G_a + G_c + G_g + G_u) + G_aa.sum()
    # Stoichiometric coefficients
    # Nucleotides [mmol.ntps/g.virus]
    S_atp = 1000 * (V_a/M_v)
    S_ctp = 1000 * (V_c/M_v)
    S_gtp = 1000 * (V_g/M_v)
    S_utp = 1000 * (V_u/M_v)
    # Amino acids [mmol.aa/g.virus]
    S_aa    = np.zeros((20,1))
    for ii in range(len(aaMets)):
        S_aa[ii] = 1000 * (V_aa[ii]/M_v)
    # Energy requirements
    # Genome: Phosphodiester bond formation products [Pyrophosphate]
    genTemp = (((countA + countC + countG + countU) * k_ppi) - k_ppi)
    genRep  = (((antiA + antiC + antiG + antiU) * k_ppi) - k_ppi)
    genTot  = genTemp + genRep
    V_ppi   = genTot
    S_ppi   = 1000 * (V_ppi/M_v)
    # Protome: Peptide bond formation [ATP + H2O]
    # Note: ATP used in this process is denoated as ATPe/Ae [e = energy version]
    spAe    = ((structCount.sum() * k_atp) - k_atp)
    npAe    = ((nonstructCount.sum() * k_atp) - k_atp)
    ppTot   = (Csp * spAe) + (Cnp * npAe)
    V_Ae    = ppTot
    S_Ae    = 1000 * (V_Ae/M_v)

    # [5] VBOF Reaction formatting and output
    # Left-hand terms: Nucleotides
    # Note: ATP term is a summation of genome and energy requirements
    S_ATP   = (S_atp + S_Ae) * -1
    S_CTP   = S_ctp * -1
    S_GTP   = S_gtp * -1
    S_UTP   = S_utp * -1
    # Left-hand terms: Amino Acids
    S_AA    = S_aa * -1
    S_AAf   = dict()
    for ii in range(len(aaMets)):
        S_AAf[aaMets[ii]] = S_AA[ii,0]
    # Left-hand terms: Energy Requirements
    S_H2O   = S_Ae * -1
    # Right-hand terms: Energy Requirements
    S_ADP   = S_Ae
    S_Pi    = S_Ae
    S_H     = S_Ae
    S_PPi   = S_ppi
    # Create reaction output
    reaction_name       = virusName + '_prodrxn_VN'
    virus_reaction      = Reaction(reaction_name)
    virus_reaction.name = virusFull + ' production reaction (created: ViraNet(c))'
    virus_reaction.subsystem                = 'Virus Production'
    virus_reaction.lower_bound              = 0
    virus_reaction.upper_bound              = 1000
    virus_reaction.objective_coefficient    = 0
    virus_reaction.add_metabolites(({
        atp_c: S_ATP,
        ctp_c: S_CTP,
        gtp_c: S_GTP,
        utp_c: S_UTP,
        ala_c: S_AAf['A'],
        arg_c: S_AAf['R'],
        asn_c: S_AAf['N'],
        asp_c: S_AAf['D'],
        cys_c: S_AAf['C'],
        gln_c: S_AAf['Q'],
        glu_c: S_AAf['E'],
        gly_c: S_AAf['G'],
        his_c: S_AAf['H'],
        ile_c: S_AAf['I'],
        leu_c: S_AAf['L'],
        lys_c: S_AAf['K'],
        met_c: S_AAf['M'],
        phe_c: S_AAf['F'],
        pro_c: S_AAf['P'],
        ser_c: S_AAf['S'],
        thr_c: S_AAf['T'],
        trp_c: S_AAf['W'],
        tyr_c: S_AAf['Y'],
        val_c: S_AAf['V'],
        h2o_c: S_H2O,
        adp_c: S_ADP,
        pi_c:  S_Pi,
        h_c:   S_H,
        ppi_c: S_PPi}))

    # DEBUG COMMENT
    # print("")
    # print("VBOF Reaction Information")
    # print(virus_reaction)
    # for ii in virus_reaction.metabolites:
    #     print("%9s : %s" % (ii.id, ii.values))
    # print(virus_reaction)

    # [6] Output return variables
    return (virus_reaction)

###############################################################################

###############################################################################
# Generation of the human-virus intergrated model (HVM): generated from a
# user-supplied host stoichiometric model and a given VBOF (generated: VBOF.py)
# Inputs:
# Model         User-supplied model file for desired host   [.mat,.xml]
# VBOF          Virus biomass objective function

# Outputs:
# hvm           Integrated host-virus model

def HVM(Model,VBOF):
    "Generate_HVM"

    # [1] Intial Setup
    # File-type dependent loading
    if ".mat".lower() in Model.lower():
        hostModel = cobra.io.load_matlab_model(Model)
    elif ".xml".lower() in Model.lower():
        hostModel  = cobra.io.read_sbml_model(Model)
    else:
        raise ValueError('Unsupported file type, unable to load model: see README')
    # [2] VBOF integration
    hostModel.add_reaction(VBOF)
    # [3] Alter the bounds of VBOF
    hostModel.reactions[-1].reversibility = False
    return (hostModel)
