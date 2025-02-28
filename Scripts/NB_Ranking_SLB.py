#!/usr/bin/env python3
# Neuroblastoma Gene Ranking script based on the NRC Dataset (GSE85047)
# Author: Christophe Van Neste and adapted by Sarah-Lee Bekaert

import pickle, re, sys, os, operator
import pandas as pd, numpy as np, matplotlib
from collections import OrderedDict
from bidali.seqanalysis import literatureLinkSearch 
from bidali.survivalinks import geneImpactSurvival
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr,pearsonr,binomtest
from leopard import Report,Section
from io import StringIO
from bidali.fegnome import enrichometer, rankSumProbability, RankSumResult
from bidali.visualizations import drawCNAcircos, dosageViolin

# Loading resources 
from lostdata.dealer.ensembl import get_ensemblGeneannot
from bidali.seqanalysis import loadHumanGenome, get_centromeres
from lostdata.dealer.entrez import get_msigdb6

# Loading data
import lostdata as LSD
from lostdata.formats import Dataset
from lostdata.dealer.home.cnvdata import get_UHRprofiles    # Fixed by adding a new import line lostdata.processing
from lostdata.dealer.home.models import get_THMYCN  # Fixed based on ConfigParser's .get() method

#mdb = get_msigdb6() # for gene set enrichtment - not needed here 
centromereshg38 = get_centromeres() #Dataset not locally available, generating: 13min
# Loads pre-computed linear model coefficients -> done by Anneleen Beckers? 
thcoeffs = get_THMYCN().lmCoeffs # 2015_Hyperplasia_ABC gedownload van Dropbox -> nieuwe input folder gemaakt and hard coded in models.py script

# Input file from Output file from expression_lm.py script
dlms = pd.read_csv(
    os.path.expanduser(
        '/code/nb_ranking/OutputData/gene_resuls_lm_SLB_V4_hg38.csv'
    ),index_col=0)

# Load CNV profiles -> frequency plot
cnvd = get_UHRprofiles() #cnvd wordt gebruikt in analyzeCNAspread en rankNBaccordingToMYCNstatus 
# FIXED! 

report = Report(
    title="Results ranking CNV genes in neuroblastoma high risk patients",
    intro= '',
    conclusion='',
    outfile=os.path.expanduser('/code/nb_ranking/OutputData')
)

report.append('Profiles overview',
              '''Amplification annotated regions are considered as regular gains
Stages: 1->1, 2->2a, 3->2b, 4->3, 5->4, 6->4S.
Stage 1 and 6 presented with MYCNamp.''')

# Creates a table showing the count of CNV profiles grouped by annotation types
report.lastSection.tabs['Annotation sizes'] = pd.DataFrame(
    cnvd.profiles.groupby('annotation').size(),columns=['size']
)
# Converts annotations labeled "amplification" to "gain." -> dus alle amplification zijn ook gewoon gains 
cnvd.profiles.annotation = cnvd.profiles.annotation.apply(
    lambda x: 'gain' if x == 'amplification' else x
)
# Updates the annotation size table after converting "amplification" to "gain."
report.lastSection.tabs['Adjusted sizes'] = pd.DataFrame(
    cnvd.profiles.groupby('annotation').size(),columns=['size']
)
# Adds a table showing the number of samples for each clinical stage.
report.lastSection.tabs['Stage distribution'] = pd.DataFrame(
    cnvd.metadata.groupby('Stage').size()
)
#cnvd.profiles = cnvd.profiles[cnvd.profiles.annotation != 'amplification']
# Calculates the proportion of CNV profiles for each chromosome
report.lastSection.tabs['Profile chromosomal distribution'] = (
    cnvd.profiles.groupby('chromosome').size()/
    len(cnvd.profiles)).sort_values(ascending=False
)

# Copy Number Alterations to determine the frequency -> figuur uit Review Bieke? Op basis van PDP! 
def analyzeCNAspread(data,mycnstatus=None):
    """
    Analyzes the penetrance of gained/lost arms across patients
    """
    import operator
    # DataFrame with CNA profiles from PDP
    profiles = data.profiles
    if mycnstatus is not None:
        samples = data.metadata
        mycnselectionset = set(samples[samples.MYCN == mycnstatus].Name)
        profiles = profiles[profiles.profile_id.isin(mycnselectionset)]
    # Exclude Ambiguous Chromosomal Arms
    profiles = profiles[~profiles.chrarm.apply(operator.contains,args=['p+q'])]
    # Unique Chromosomal Arms Per Patient
    patientarms = profiles[['profile_id','chrarm']]
    patientarms = patientarms.drop_duplicates()
    # Groups data by chromosomal arm and calculates the proportion of patients affected for each arm.
    chrarmpenetrance = (patientarms.groupby('chrarm').size()/
                        len(set(patientarms.profile_id))).sort_values()
    return chrarmpenetrance # Output bv 17q - 80%

# Ranks chromosomal arms and genes in neuroblastoma patients based on their MYCN amplification status.
# Space make the function break up - python is very specific about his tabs/spaces
def rankNBaccordingToMYCNstatus(mycnstatus,data,minPenetrance=.25,penetranceOnly=True,
                                survivalWithinGroupOnly=False,filterNormalGainLosses=False,
                                filterOnGenePenetrance=False):
    """
    mycnstatus => 0 == non amplified, 1 == amplified
    min[imal]Penetrance, has been set to 0.25, within the PDP dataset this includes all common
    non MYCN amplified gained/lost chromosome arms. This is also applied to the gene level:
    genes that are only present in less than minPenetrance of its chr arm gain/losses are discarded
    penetranceOnly: implies segments are combined from the same sample, ergo a gene can count only
    once for a sample to determine its percentage
    """
    #Setup report section
    section = Section(
        title='Ranking MYCN {} cases'.format('amp' if mycnstatus else 'sc'),
        text='minPenetrance={},survivalWithinGroupOnly={}\n '.format(minPenetrance,survivalWithinGroupOnly),
        clearpage=True)
    profiles = data.profiles
    samples = data.metadata
    #Chromosomes to rank
    chralterationPercs = analyzeCNAspread(data,mycnstatus).sort_values(ascending=False)
    section.tabs['Chromosome arm alteration sample spread'] = chralterationPercs
    selectedAlterations = chralterationPercs.index[chralterationPercs>minPenetrance]
    #chralterationPercs = chralterationPercs.filter(regex='[^+][pq]$')
    mycnselectionset = set(samples[samples.MYCN == mycnstatus].Name)
    profiles = profiles[profiles.profile_id.isin(mycnselectionset)]
    ## Chr arm alterations versus chr arm global characterisitcs
    centromereshg38['plen'] = centromereshg38.left_base
    df = pd.DataFrame(chralterationPercs,columns=['chralts'])
    df['genes'] = df.apply(lambda x: centromereshg38.loc[x.name[:-1]][x.name[-1]+'_genes'],axis=1)
    df['len'] = df.apply(lambda x: centromereshg38.loc[x.name[:-1]][x.name[-1]+'len'],axis=1)
    df['gdensity'] = df.genes/df.len
    ## Gains over losses ratio
    gols = profiles.groupby('chrarm').apply(lambda x: x.groupby('annotation').size())
    if type(gols) == pd.Series: gols = gols.unstack()
    gols['ratio'] = gols.gain/gols.loss
    gols['binomtest'] = gols['binomtest'] = gols.fillna(0).T.apply(lambda x: binomtest(int(x.gain),int(x.gain+x.loss)).pvalue)
    #gols['binomtest'] = gols.fillna(0).T.apply(lambda x: binomtest(int(x.gain),int(x.gain+x.loss)).pvalue)
    unsignols = gols.dropna()[gols.dropna().binomtest > 0.05].index #chr arms without either losses or gains should not be included
    lossarms = gols[gols.ratio<1].index
    # Retrieve genes
    genannot = get_ensemblGeneannot()
    ## Process each important chr arm separately
    chresults = OrderedDict()
    for chrarm in selectedAlterations:
        chraprofs = profiles[profiles.chrarm==chrarm]
        togain = chraprofs.groupby('annotation').size()
        print(togain)
        try: togain = togain.loc['gain'] > togain.loc['loss']
        except KeyError: togain = 'gain' in togain.index
        chraprofs = chraprofs[chraprofs.annotation == ('gain' if togain else 'loss')]
        chraprofs['genes'] = chraprofs.T.apply(lambda x: {f.attributes['gene_name'][0] for f in genannot.region('{}:{}-{}'
                                            .format(x.chromosome,int(x.min38),int(x.max38)),featuretype='gene')})
        chraprofs['nrGenes'] = chraprofs.genes.apply(len)  
        if penetranceOnly:
            combinedSegments = chraprofs.groupby('profile_id').apply(lambda x: set.union(*list(x['genes'])))
            combinedmin = chraprofs.groupby('profile_id').apply(lambda x: x.min38.min())
            combinedmax = chraprofs.groupby('profile_id').apply(lambda x: x.max38.max())
            chraprofs = chraprofs[~chraprofs.profile_id.duplicated()].copy()
            chraprofs.index = chraprofs.profile_id
            chraprofs.genes = combinedSegments
            chraprofs.min38 = combinedmin
            chraprofs.max38 = combinedmax
            chraprofs.nrGenes = chraprofs.genes.apply(len)     
        chragenes = {}
        for gs in chraprofs.genes:
            for g in gs:
                try: chragenes[g]+=1
                except KeyError: chragenes[g]=1
        chragenes = pd.DataFrame({'number':chragenes})
        chragenes['percentage'] = chragenes.number/len(chraprofs)
        chragenes['penetrance'] = chragenes.percentage*chralterationPercs.loc[chrarm]
        if filterOnGenePenetrance: chragenes = chragenes[chragenes.penetrance >= minPenetrance]
        # Copy Number Frequency 
        chragenes['CNrank'] = chragenes.number.rank()
        chresults[chrarm] = ('gain' if togain else 'loss',chraprofs,chragenes)
    chrankedgenes = OrderedDict()
    for chrarm in chresults:
        togain,chraprofs,chragenes = chresults[chrarm]
        if filterNormalGainLosses:
            chrarmLength = centromereshg38.loc[chrarm[:-1]][chrarm[-1]+'len']
            normalGainLosses = fitfunctions['len'][2](chrarmLength)
            # Calculate fitted normal percentage of gain losses + 2SD
            TT=np.array(((chrarmLength,1),))
            yi = np.dot(TT, fitfunctions['len'][0])
            C_yi = np.dot(TT, np.dot(fitfunctions['len'][1], TT.T))
            sig_yi = np.sqrt(np.diag(C_yi))
            acceptedGainLosses = float(yi+2*sig_yi)
            normalChrPercentage = (acceptedGainLosses*len(profiles)/2)/len(chraprofs)
            print(chrarm,len(chragenes),sum(chragenes.percentage < normalChrPercentage),sum(chragenes.percentage >= normalChrPercentage))
            chragenes = chragenes[chragenes.percentage >= normalChrPercentage].copy()
            if len(chragenes) == 0:
                print('All',chrarm,'genes filtered according to normalChrPercentage')
                continue      
        chragenes['CNrank'] = chragenes.number.rank()
        chragenes['THMYCNdiff'] = chragenes.apply(lambda x: thcoeffs.loc[x.name.capitalize()]['lineardiff']
                                                  if x.name.capitalize() in thcoeffs.index else np.nan,axis=1)
        # lineair model difference for TH-MYCN
        chragenes['THMYCNrank'] = chragenes.THMYCNdiff.rank(ascending=togain == 'gain')
        #LM expression model -> from the NRC dataset computed model dlms
        chragenes['dosagelm'] = chragenes.T.apply(lambda x: np.nan if x.name not in dlms.index else dlms.loc[x.name]['cnfocus'+togain])
        # Dosage Sensitivty -> from the NRC dataset computed model dlms
        chragenes['dosagerank'] = chragenes.dosagelm.rank(ascending=togain == 'gain')
        # Used to rank genes based on their survival association -> from the NRC dataset computed model dlms
        chragenes['risklm'] = chragenes.T.apply(lambda x: np.nan if x.name not in dlms.index else dlms.loc[x.name]['metadatahighstage_'+('amp' if mycnstatus else 'sc')])
        # Risk association based on chr arm -> riskrank ranks genes based on their risklm values, either in ascending or descending order, depending on the chromosomal arm's alteration type (gain or loss).
        chragenes['riskrank'] = chragenes.risklm.rank(ascending=togain == 'gain')
        #Filter genes that lack certain characteristics
        print(chrarm,'before filtering lacking more than 1 characteristic',chragenes.shape[0],'genes')
        chragenes = chragenes[chragenes.isnull().sum(axis=1) <= 2]
        print('after filtering',chragenes.shape[0])
        #Results
        chrankedgenes[chrarm+'_'+togain] = chragenes
    return section,chresults,chrankedgenes, Dataset(
        minPenetrance=minPenetrance,
        df=df,
        lossarms=lossarms,
        unsignols=unsignols)

# Calculating within HR subgroups
# Processing with NRC dosage lm
# Filtering genes lacking more than 1 characteristic
# MYCNsc
section,chresultsMYCNsc,chrankedgenesMYCNsc,sc_cn_dataset = rankNBaccordingToMYCNstatus(0,cnvd,survivalWithinGroupOnly=True)
section._reportSection = report #TODO -> not proper way, when appending should be set CVN opmerking
section._parentSection = report
report.sections.append(section)

# Reranking output -> added some lines to make the dir
output_dir = os.path.expanduser("/code/nb_ranking/OutputData")
os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist
for chrarm in chrankedgenesMYCNsc:
    chrarmranking = chrankedgenesMYCNsc[chrarm][['penetrance','dosagelm','risklm']].copy()
    chrarmranking['prod_pct_rank'] = chrarmranking.rank(pct=True).prod(axis=1)
    chrarmranking.sort_values('prod_pct_rank',inplace=True,ascending=False)
    chrarmranking.to_csv(os.path.expanduser(f"/code/nb_ranking/OutputData/sc{chrarm}_ranking.csv"))


# MYCNamp                                                                                  # 1 = MYCN amplified, ... zelfde      
section,chresultsMYCNamp,chrankedgenesMYCNamp,amp_cn_dataset = rankNBaccordingToMYCNstatus(1,cnvd,survivalWithinGroupOnly=True)
section._reportSection = report
section._parentSection = report
report.sections.append(section)
 
# Ranking per chr arm alone for MYCN amplified cases
for chrarm in chrankedgenesMYCNamp:
    # Selection of the 3 parameters where ranking is based on for that specific chr arm
    chrarmranking = chrankedgenesMYCNamp[chrarm][['penetrance','dosagelm','risklm']].copy()
    # ranks each value in the selected columns relative to others, as a percentage (values between 0 and 1)
    # Hier doen we eindelijk de combined ranking! 
    chrarmranking['prod_pct_rank'] = chrarmranking.rank(pct=True).prod(axis=1)
    # Sorting the results and save 
    chrarmranking.sort_values('prod_pct_rank',inplace=True,ascending=False)
    chrarmranking.to_csv(os.path.expanduser(f"/code/nb_ranking/OutputData/amp{chrarm}_ranking.csv"))
