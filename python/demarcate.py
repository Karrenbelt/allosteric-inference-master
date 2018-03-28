import settings
import os, re, cobra
import pandas as pd
import numpy as np

### read models, core and genome scale (iJO1366) 
GSMM = cobra.io.read_sbml_model(settings.ECOLI_SBML_FNAME)
core = cobra.io.read_sbml_model(settings.ECOLI_CORE_FNAME)

## 1. convert core model to iJO1366 name space, mark descrepancies between models kept

# remove boundary species, don't occur in iJO1366
for rxn in core.reactions:
    for met in rxn.metabolites:
        if met.id.endswith('_b'):
            if rxn.reversibility == False:
                rxn.reaction = re.split(' --> ',rxn.reaction)[0]+' --> '
            else:
                rxn.reaction = re.split(' <=> ',rxn.reaction)[0]+' <=> '
            core.metabolites.remove(met) # does not remove core mets? 
            
# change inconsistence in metabolite naming, adopt iJO1366 standard
manual_rxn_map = {'EX_glc_e':'EX_glc__D_e','EX_gln_L_e':'EX_gln__L_e','EX_glu_L_e':'EX_glu__L_e','EX_lac_D_e':'EX_lac__D_e','EX_mal_L_e':'EX_mal__L_e'}
for met in core.metabolites:
    if met not in GSMM.metabolites:
        met.id = met.id.replace('_','__',1)
for rxn in core.reactions:
    if rxn not in GSMM.reactions:
        if rxn.id in manual_rxn_map:
            rxn.id = manual_rxn_map[rxn.id]
core.repair()


## 2. add additional reactions of interest from iJO1366

# # add transport reactions in similar style as the others in the core model
# def add_from_iJO1366(rxn_ids_list):
#     [core.add_reactions([GSMM.reactions.get_by_id(rxn)]) for rxn in rxn_ids_list]
#     return core
# ### maybe interesting but not added:
# # + ['FBA3', 'PFK_3'] + ['MGSA','LGTHL','GLYOX'] + ['OAADC']
# # [CITL','FRD2','FRD3','MDH2','MDH3'] + ['F6PA']

# ## adding pathways for gluconate
# GLCNt2r = cobra.Reaction('GLCNt2r') # GLCNt2rpp, GLCNtex
# h_e = GSMM.metabolites.get_by_id('h_e')
# h_c = GSMM.metabolites.get_by_id('h_c')
# glcn_e = GSMM.metabolites.get_by_id('glcn_e')
# glcn_c = GSMM.metabolites.get_by_id('glcn_c')
# GLCNt2r.add_metabolites({glcn_e: -1.0, h_e:-1.0, glcn_c: 1.0, h_c:1.0})
# core.add_reactions([GLCNt2r])
# core = add_from_iJO1366(['GNK','EX_glcn_e'])

# ## adding pathways for glycerol
# GLYCt = cobra.Reaction('GLYCt') # GLYCtex, GLYCtpp
# glyc_e = GSMM.metabolites.get_by_id('glyc_e')
# glyc_c = GSMM.metabolites.get_by_id('glyc_c')
# GLYCt.add_metabolites({glyc_e: -1.0, glyc_c: 1.0})
# core.add_reactions([GLYCt])
# core = add_from_iJO1366(['GLYCDx','DHAPT','G3PD2','G3PT','GLYK','EX_glyc_e'])

# ## adding pathways for galactose
# GALabc = cobra.Reaction('GALabc') # GALtex, GALabspp
# gal_e = GSMM.metabolites.get_by_id('gal_e')
# gal_c = GSMM.metabolites.get_by_id('gal_c')
# atp_c = GSMM.metabolites.get_by_id('atp_c')
# adp_c = GSMM.metabolites.get_by_id('adp_c')
# GALabc.add_metabolites({gal_e: -1.0, atp_c: -1.0, gal_c: 1.0, adp_c: 1.0})
# core.add_reactions([GALabc])

# UGLT = cobra.Reaction('UGLT') # effectively removes udpg_c and updgal_c from the original reaction
# gal1p_c = GSMM.metabolites.get_by_id('gal1p_c') 
# g1p_c = GSMM.metabolites.get_by_id('g1p_c')
# UGLT.add_metabolites({gal1p_c: -1.0, g1p_c: 1.0})
# core.add_reactions([UGLT])
# core = add_from_iJO1366(['PGMT','GALKr','EX_gal_e'])

# ## changing pathways for pyruvate
# core.remove_reactions(['FRUpts2']) # change entry point to fpb_c instead of f6p_c (ambiguous), to match Gerosa 2015
# FRUpts = cobra.Reaction('FRUpts')
# fru_e = GSMM.metabolites.get_by_id('fru_e')
# f1p_c = GSMM.metabolites.get_by_id('f1p_c')
# FRUpts.add_metabolites({fru_e: -1.0, f1p_c: 1.0})
# core.add_reactions([FRUpts])
# core = add_from_iJO1366(['FRUK'])

# ## add ED pathway, ACS for acetate, and MOX for malate oxidase (MQO, MDH2, MDH3: all the same reaction, diff cofactor)
# core = add_from_iJO1366(['EDD','EDA'] + ['ACS'] + ['MOX'])

# core.reactions.get_by_id('EX_glc__D_e').bounds = [0, 1000] # is not set to zero in the core model initially

# core.repair()

# cobra.io.write_sbml_model(core, settings.CACHE_DIR+'/extended_core.xml')



### differences kept:
# the GSMM contains a periplasmic compartment, core model does not.
# ETOHtrpp + ETOHtex ~= ETOHt2r (proton symport in core, not in iJO1366)
# FRD2, FRD3 ~= FRD7 (mql8_c, 2dmmql8_c VS q8h2_c)
# BIOMASS_Ec_iJO1366_WT_53p95M ~= Biomass_Ecoli_core_w_GAM
# approximate_map = {'ETOHt2r':'ETOHtrpp','FORti':'FORtppi','FRD7':'FRD2','Biomass_Ecoli_core_w_GAM':'Biomass_Ecoli_core_w_GAM',
# additionally all transport and exchange reactions in the core model have been merged (e.g. SUCCt3pp + SUCCtex --> SUCCt3)

# missing in Glycolysis core model:
# F6PA   f6p_c <=> dha_c + g3p_c
# G1PPpp     g1p_p + h2o_p --> glc__D_p + pi_p
# G6PP   g6p_c + h2o_c --> glc__D_c + pi_c
# GLBRAN2    glycogen_c --> bglycogen_c
# GLCP   glycogen_c + pi_c --> g1p_c
# GLCP2      bglycogen_c + pi_c --> g1p_c
# GLCS1      adpglc_c --> adp_c + glycogen_c + h_c
# GLDBRAN2   bglycogen_c --> glycogen_c
# GLGC   atp_c + g1p_c + h_c --> adpglc_c + ppi_c
# HEX1   atp_c + glc__D_c --> adp_c + g6p_c + h_c

# missing in PPP core model:
# EDA    2ddg6p_c --> g3p_c + pyr_c
# EDD    6pgc_c --> 2ddg6p_c + h2o_c
# FBA3   s17bp_c <=> dhap_c + e4p_c
# PFK_3      atp_c + s7p_c --> adp_c + h_c + s17bp_c

# missing in TCA core model:
# CITL   cit_c --> ac_c + oaa_c
# FRD2   fum_c + mql8_c --> mqn8_c + succ_c
# FRD3   2dmmql8_c + fum_c --> 2dmmq8_c + succ_c
# MDH2   mal__L_c + q8_c --> oaa_c + q8h2_c
# MDH3   mal__L_c + mqn8_c --> mql8_c + oaa_c
# MOX    mal__L_c + o2_c <=> h2o2_c + oaa_c

# missing in Anaplerosis core model:
# PPA    h2o_c + ppi_c --> h_c + 2.0 pi_c
# PPA2   h2o_c + pppi_c --> h_c + pi_c + ppi_c
