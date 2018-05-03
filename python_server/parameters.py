#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 11:43:03 2018

@author: Zarathustra
"""

import os, settings, cobra, pickle, json, re, boolean
import bigg, kegg
import numpy as np
import pandas as pd
from cobra.util import create_stoichiometric_matrix
from scipy import constants
from uncertainties import ufloat

condition = 'Succinate'
#%% Parsing data from Link

df = pd.read_excel(settings.DATA_DIR+'/Link_2014_S2.xlsx')
df.columns = list(range(0,len(df.columns)))
df = df[3:-11]
met_names = [df.iloc[i][1] for i in range(0,len(df.index),2)]
rel_abs = df[0][::2]

pyr_mean  = [df.iloc[i][2:9].values for i in range(0,len(df.index),2)]
pyr_sd    = [df.iloc[i+1][2:9].values for i in range(0,len(df.index),2)]
fru_mean  = [df.iloc[i][9:-1].values for i in range(0,len(df.index),2)]
fru_sd    = [df.iloc[i+1][9:-1].values for i in range(0,len(df.index),2)]
ctrl_mean = [df.iloc[i][-1:].values for i in range(0,len(df.index),2)]
ctrl_sd   = [df.iloc[i+1][-1:].values for i in range(0,len(df.index),2)]

# df_pyr_mean = pd.DataFrame(pyr_mean, index=met_names, columns=time)
columns = [10,15,25,40,45,55,70]
df_pyr_mean = pd.DataFrame(pyr_mean, columns=columns, index = met_names)
df_pyr_sd   = pd.DataFrame(pyr_sd,   columns=columns, index = met_names)
df_fru_mean = pd.DataFrame(fru_mean, columns=columns, index = met_names)
df_fru_sd   = pd.DataFrame(fru_sd,   columns=columns, index = met_names)

mapping = {    
        'asparagine':'asn__L_c',
        'glutamine':'gln__L_c',
        'citrulline':'citr__L_c',
        'diaminopimelate':'diaminopimelate',
        'Homoserine':'hom__L_c',
        'Ga6P':'gam6p_c',
        'Guanine':'gua_c',
        'Tyrosine':'tyr__L_c',
        'phenylalanine':'phe__L_c', 
        'aspartate':'asp__L_c',
        'glutamate':'glu__L_c',
        'Lactate':'lac__D_c', 
        'tryptophane':'trp__L_c',
        'G6P':'g6p_c',
        'Ms6P':'man6p',
        'R5P':'r5p_c',
        'F6P':'f6p_c',
        'GlycerolP':'glyc1p_c',
        'G1P':'g1p_c',
        'S7P':'s7p_c',
        'Ru5P':'ru5p__D_c',
        'Xu5P':'xu5p__D_c', 
        'GTTred':'gthrd_c',
        'DHAP':'dhap_c',
        'NAD':'nad_c',
        'Panthothenate':'pnto__R_c',
        'cAMP':'camp_c',
        'Succinate':'succ_c',
        'GTTox':'gthox_c',
        'Malate':'mal__L_c',
        'UDP-hexose':'UDP-hexose',
        'alpha ketoglutarate':'akg_c',
        'fumarate':'fum_c',
        'ADP-hexose':'ADP-hexose', 
        'ADP-pentose':'ADP-pentose',
        'xPG':'3pg_c',
        '6PG':'6pgc_c', 
        'NADH':'nadh_c',
        'Aconitate':'acon_C_c',
        'NADP':'nadp_c',
        'PEP':'pep_c',
        'FBP':'fdp_c', 
        'methylcitrate_methylisocitrate':'2mcit_c',
        'citrate_isocitrate':'cit_c', ### citrate = citrate_isocitrate - isocitrate
        'isocitrate':'icit_c',
        'FMN':'fmn_c',
        'BPG':'13dpg_c',
        'FAD':'fad_c',
        'AcCoA':'accoa_c',
        'NADPH':'nadph_c',
        'PRPP':'prpp_c',
        'ATP':'atp_c',
        'ADP':'adp_c',
        'AMP':'amp_c',
        'GTP':'gtp_c',
        'GDP':'gdp_c',
        'GMP':'gmp_c',
}

df_pyr_mean.loc['citrate_isocitrate'] = df_pyr_mean.loc['citrate_isocitrate'] - df_pyr_mean.loc['isocitrate']
df_fru_mean.loc['citrate_isocitrate'] = df_fru_mean.loc['citrate_isocitrate'] - df_fru_mean.loc['isocitrate']
df_pyr_mean.rename(index=mapping)
df_pyr_sd.rename(index=mapping)
df_fru_mean.rename(index=mapping)
df_fru_sd.rename(index=mapping)

df_pyr_mean.to_csv(settings.CACHE_DIR+'/ts_pyr_mean_link2013.csv')
df_pyr_sd.to_csv(settings.CACHE_DIR+'/ts_pyr_sd_link2013.csv')
df_fru_mean.to_csv(settings.CACHE_DIR+'/ts_fru_mean_link2013.csv')
df_fru_sd.to_csv(settings.CACHE_DIR+'/ts_fru_sd_link2013.csv')

# df1 = pd.read_csv(settings.CHACHE_DIR+'/ecoli_metabolites_link2013.csv', index_col=0)
df2 = pd.read_csv(settings.DATA_DIR+'/ecoli_metabolites_kochanowski2017.csv', index_col=0)
df3 = pd.read_csv(settings.DATA_DIR+'/ecoli_metabolites_gerosa2015.csv', index_col=0)

#%% Make DataFrame Template

def string_split(s, c, n):
    """
        Uses regex to split a string on the n-th occurence of a character
        args:
            s (string) - string to be split
            c (string) - character(s) on which to split
            n (int)    - n-th occurence
        returns:
            l (tuple)  - the split string
    """
    regex = r'^((?:[^%c]*%c){%d}[^%c]*)%c(.*)' % (c,c,n-1,c,c)
    m = re.match(regex, s)
    if m: return m.groups()
    else: return ()

def get_rxn_and_met_id(string_name):
    rxn_id = re.findall(r'R[0-9]+', string_name)
    met_id = re.findall(r'S[0-9]+', string_name)
    if not rxn_id: rxn_id = [np.nan]
    if not met_id: met_id = [np.nan]
    assert(len(rxn_id)==len(met_id)==1)
    return rxn_id[0], met_id[0]

def pb_data_template(species_and_parameters):
    """
        Generates a template for parameter balancing
        args:
            species_and_parameters (dictionary) - model species and parameters
        returns:
            template (DataFrame) - The fields 
            !Reaction, !Compound, !Unit, '!QuantityType' are already filled out
            !Mean, !Std, !Min, !Max remain empty
    """
    header  = '!!SBtab parameters extended core model E. coli'
    columns = ['!QuantityType','!Reaction','!Compound','!Mean','!Std','!Unit','!Min','!Max']
    
    index = []
    for key in species_and_parameters.keys():
        if key.count('_') == 3: # compound Vmax ID (e.g. E_R0_kcat_R0)
            index.extend(list(string_split(key,'_',2)))
        else:
            index.append(key)
        
    df = pd.DataFrame(index=index, columns=columns)
    for item in df.index:
        rxn_id, met_id = get_rxn_and_met_id(item)
        df.loc[item,'!Reaction']         = rxn_id # mapping.get(rxn_id)
        df.loc[item,'!Compound']         = met_id # mapping.get(met_id)
        if item.startswith('kcat_'):
            df.loc[item,'!QuantityType'] = 'catalytic rate constant'
            df.loc[item,'!Unit']         = '1/s'
        elif item.startswith('E_'):
            df.loc[item,'!QuantityType'] = 'concentration of enzyme'
            df.loc[item,'!Unit']         = 'molecules/cell'
        elif item.startswith('Km_'):
            df.loc[item,'!QuantityType'] = 'Michaelis constant'
            df.loc[item,'!Unit']         = 'mM'
        elif item.startswith('Keq_'):
            df.loc[item,'!QuantityType'] = 'equilibrium constant'
            df.loc[item,'!Unit']         = 'dimensionless'
        elif item.startswith('Ki_'):
            df.loc[item,'!QuantityType'] = 'inhibitory constant'
            df.loc[item,'!Unit']         = 'mM' 
        elif item.startswith('v_'):
            df.loc[item,'!QuantityType'] = 'flux'
            df.loc[item,'!Unit']         = 'mmol/gCDW/h'

    return df


model_mapping = json.load(open(settings.CACHE_DIR+'/mapping.json'))
species_and_parameters = json.load(open(settings.CACHE_DIR+'/model_species_and_parameters.json'))
df = pb_data_template(species_and_parameters)

#%% Parsing enzyme data
         
def recursive(boolean_rules, enzyme_cond_data):
    
    def aggregate(agg, boolean_rule_list):
        values = []
        for boolean_rule in boolean_rule_list:
            value = recursive(boolean_rule, enzyme_cond_data)
            if value is not None:
                values.append(value)
        if values != []:
            return agg(values)
        else:
            return None

    if type(boolean_rules) == boolean.boolean.Symbol: 
        locus_tag = str(boolean_rules)
        if locus_tag in enzyme_cond_data.index:
            value = enzyme_cond_data[locus_tag]
            if type(value) is not np.int64: # multiple values for the same locus_tag
                pass # this need manual curation: AND or OR logic applies?
            else:
                return value
        else: 
            return None
        
    elif type(boolean_rules) == boolean.boolean.OR:
        return aggregate(sum, boolean_rules.args)
    elif type(boolean_rules) == boolean.boolean.AND:
        return aggregate(min, boolean_rules.args)
    else:
        raise ValueError('exception: type {}'.format(type(boolean_rules)))
            
    
enzyme_data = pd.read_csv(settings.ECOLI_PROT_FNAME).set_index('Bnumber')
GSMM = cobra.io.read_sbml_model(settings.ECOLI_GSMM_FNAME)   
extended_core = cobra.io.read_sbml_model(settings.ECOLI_EXCORE_FNAME)
algebra = boolean.BooleanAlgebra()

enzyme_cond_mean = enzyme_data[condition+' (mean)']
#enzyme_cond_std  = enzyme_data[condition+' (cv)'].multiply(enzyme_cond_mean)/100

enzyme_cond_data = enzyme_cond_mean
for rxn in extended_core.reactions:
    gene_rules = rxn.gene_reaction_rule
    enzyme_count = []
    if len(gene_rules) != 0:
        boolean_rules = algebra.parse(gene_rules)
        enzyme_count = recursive(boolean_rules, enzyme_cond_data)
    else: 
        enzyme_count = np.nan
        
    if not rxn.id.startswith('EX_'):
        df.loc['E_'+model_mapping.get(rxn.id),'!Mean'] = enzyme_count
        


#%% Parsing Michaelis-Menten constants

#creates list with bigg.reactionIDs corresponding to the EC_numbers
def mapping_EC_to_bigg(dataframe):
    bigg = []
    for met in list(dataframe['EC_number'].values):
        if mapping['bigg.reaction'].loc[mapping['EC_number']==met].values.size != 0:
            bigg.append(str(mapping['bigg.reaction'].loc[mapping['EC_number']==met].values[0]))
        else:
            bigg.append(None)
    return bigg

km = pd.read_csv(settings.BRENDA_KM, index_col=0, encoding='ISO-8859-1')
value_col = 'KM_Value'
organism = 'Escherichia coli'
km_coli = km[km['Organism'].str.lower() == organism.lower()]
# print((km_coli[value_col] > 0).sum())
km_filtered = km_coli[(pd.isnull(km_coli['Commentary'])) |
      ((km_coli['Commentary'].str.find('mutant') == -1) &
       (km_coli['Commentary'].str.find('mutation') == -1) &
       (km_coli['Commentary'].str.find('variant') == -1) &
       (km_coli['Commentary'].str.find('genetically engineered') == -1))]

km = km_filtered[pd.notnull(km_filtered['bigg.metabolite'])].copy()
km['bigg.metabolite'] = km['bigg.metabolite'].str.lower()
km_val = km[km[value_col] > 0]
print(km_val.shape[0])
print(km_val.groupby(('bigg.metabolite', 'EC_number')).first().shape[0])
print(km_val.groupby('bigg.metabolite').first().shape[0])
print(km_val.groupby('EC_number').first().shape[0])

km_val.groupby(('bigg.metabolite', 'EC_number')).first()
bigg_mapping = bigg.BiGG._get_reaction_df()
bigg_ec_dict = dict(zip(bigg_mapping.loc[:,'bigg.reaction'], 
         bigg_mapping.loc[:,'EC_number']))
bigg_ec_dict.update({v: k for k, v in bigg_ec_dict.items()})


km = km.set_index('EC_number')
for idx in df.index:
    if idx.startswith('Km_'):
        rxn, met = get_rxn_and_met_id(idx)
        mapped_met = model_mapping[met][:2] # strip compartment
        ec_number = bigg_ec_dict.get(model_mapping[rxn].lower())
        if ec_number in km.index: # if there exists a ec_number, and in data
            km_vals = km.loc[ec_number,:]
            
            if isinstance(km_vals, pd.Series):
                km_vals = pd.DataFrame(km_vals).T
            
            kms_rxn = km_vals.set_index('bigg.metabolite')
            if mapped_met in kms_rxn.index: # if met in data
                values = kms_rxn.loc[mapped_met,'KM_Value']
                uniques = list(set(values))
                assert(len(uniques) == 1)
                df.loc[idx, '!Mean'] = uniques[0]
            

#%% Parsing data for catalytic constants

kcat = pd.read_csv('kcat1s.csv', sep=';') # <--- where does this come from?!
kcat.columns = kcat.iloc[1]
kcat = kcat.drop(kcat.index[0])
kcat = kcat.drop(kcat.index[0])
kcat.index = kcat['reaction (model name)']
kcat.rename(columns={'kcat per active site [1/s]':'!Mean'}, inplace=True)

for i in kcat.index:
    ID = 'kcat_' + i
    if ID in df.index:
        df.at[ID, '!Mean'] = kcat.at[i,'!Mean'] 
        
#%% Parsing data for quilibrium constants

def overflow_values(liste):
    new_liste = [10**6 if str(s) == 'inf' else s for s in liste]
    return new_liste

therm = pd.read_csv(settings.ECOLI_THERMO_FNAME, index_col = 0)
value_col = 'dG0_prime'
organism = 'Escherichia coli'

T = 300
R = 8.314472*10**(-3) # R = 8.314 J mol-1 K-1 or 0.008314 kJ mol-1 K-1.
Keq = np.exp((-therm['dG0_prime'])/(R*T))
std = Keq*therm['dG0_prime_std']/(R*T)

Keq = overflow_values(Keq)
std = overflow_values(std)

new = pd.DataFrame(index=therm.index, columns=['Keq','std'])
new['Keq'] = Keq
new['std'] = std

for item in df.index:
    if item.startswith('Keq_'):
        rxn = item.split('_',1)[1].lower()
        if rxn in new.index:
            df.loc[item,'!Mean'] = new.loc[rxn,'Keq']
            df.loc[item,'!Std'] = new.loc[rxn,'std']

therm2 = pd.read_csv('ThermoDataiJO1366.csv', index_col = 0,sep = ';')
Keq2 = np.exp((-therm2['dG_r_prime_0'].fillna(0))/(R*T))
std2 = Keq2*therm2['stderr'].fillna(0)/(R*T)

Keq2 = overflow_values(Keq2)
std2 = overflow_values(std2)

new2 = pd.DataFrame(index=therm2.index.str.lower().str.strip(), columns=['Keq','std'])
new2['Keq'] = Keq2
new2['std'] = std2

for item in df.index:
    if item.startswith('Keq_'):
        rxn = item.split('_',1)[1].lower()
        if rxn in new2.index:
            print(new2.loc[rxn,'Keq'])
#             print(df.loc[item,'!Mean'], new.loc[rxn,'Keq'])
#             df.loc[item,'!Mean'] = new.loc[rxn,'Keq']
#             df.loc[item,'!Std'] = new.loc[rxn,'std']
