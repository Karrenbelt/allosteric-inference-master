import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import imp
import json
import inspect
import urllib
from contextlib import closing

CONDITIONS = ['Glucose','Gluconate','Galactose','Pyruvate','Glycerol','Succinate','Acetate','Fructose']
COND2BIGG  = {'Glucose':'EX_glc__D_e','Pyruvate':'EX_pyr_e','Fructose':'EX_fru_e','Acetate':'EX_ac_e', 
              'Succinate':'EX_succ_e','Glycerol':'EX_glyc_e','Gluconate':'EX_glcn_e','Galactose':'EX_gal_e'}

SCRIPT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
BASE_DIR = os.path.join(*os.path.split(SCRIPT_DIR)[0:-1])

DATA_DIR = os.path.join(BASE_DIR, 'data')
CACHE_DIR = os.path.join(BASE_DIR, 'cache')
RESULT_DIR = os.path.join(BASE_DIR, 'result')
BRENDA_DIR = os.path.join(BASE_DIR, 'data/brenda_query_2016_06_20')

## ID mapping files
KEGG2CHEBI_FNAME         = os.path.join(DATA_DIR, 'kegg2chebi.csv')
BIGG_METABOLITE_FNAME    = os.path.join(DATA_DIR, 'bigg_models_metabolites.txt')
BIGG_REACTION_FNAME      = os.path.join(DATA_DIR, 'bigg_models_reactions.txt')
ECOLI_THERMO_CACHE_FNAME = os.path.join(DATA_DIR, 'ecoli_thermodynamics.csv')

## E. coli model files
ECOLI_JSON_FNAME = os.path.join(DATA_DIR, 'iJO1366.json')
ECOLI_SBML_FNAME = os.path.join(DATA_DIR, 'iJO1366.xml')
ECOLI_CORE_FNAME = os.path.join(DATA_DIR, 'Orth2011_core.xml')
ECOLI_XLS_FNAME  = os.path.join(DATA_DIR, 'inline-supplementary-material-2.xls')

## Data files
# GEROSA_FLUX_FNAME = os.path.join(DATA_DIR, 'ecoli_fluxes_gerosa2015.csv')
ECOLI_METAB_FNAME = os.path.join(DATA_DIR, 'ecoli_metabolites_gerosa2015.csv')
ECOLI_METAB_FNAME = os.path.join(DATA_DIR, 'ecoli_metabolites_kochanowski2017.csv')
ECOLI_PROT_FNAME  = os.path.join(DATA_DIR, 'ecoli_proteins_schmidt2015.csv')
ECOCYC_REG_FNAME  = os.path.join(DATA_DIR, 'EcoCycRegulation.csv')
MANUAL_REG_FNAME  = os.path.join(DATA_DIR, 'ecoli_reg_manual_curation.csv')
ECOLI_TFN_FNAME   = os.path.join(DATA_DIR, 'RegulonDB_network_tf_gene')

## Brenda raw data
RAW_BRENDA_ACTIVATING = os.path.join(BRENDA_DIR, 'activating.csv')
RAW_BRENDA_COFACTORS  = os.path.join(BRENDA_DIR, 'cofactors.csv')
RAW_BRENDA_INHIBITING = os.path.join(BRENDA_DIR, 'inhibiting.csv')
RAW_BRENDA_KI         = os.path.join(BRENDA_DIR, 'ki.csv')
RAW_BRENDA_KM         = os.path.join(BRENDA_DIR, 'km.csv')
RAW_BRENDA_METAL      = os.path.join(BRENDA_DIR, 'metal.csv')
RAW_BRENDA_TURNOVER   = os.path.join(BRENDA_DIR, 'turnover.csv')

BRENDA_KM         = os.path.join(DATA_DIR, 'km.csv')
BRENDA_TURNOVER   = os.path.join(DATA_DIR, 'turnover.csv')


## cache files
ECOLI_EXCORE_FNAME = os.path.join(CACHE_DIR, 'extended_core.xml')
ECOLI_GEROSA_FLUX  = os.path.join(CACHE_DIR, 'gerosa_fluxes_per_rxn.csv')


