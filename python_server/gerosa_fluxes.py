#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 12:00:30 2018

@author: Zarathustra
"""

import settings, cobra, re
import pandas as pd
import numpy as np
import sympy
from uncertainties import ufloat
from cobra.util import create_stoichiometric_matrix

def Physiology(filename):
    """
        Parsing of the physiology data from Gerosa (2015)
    """
    short_names = ['growth_rate','uptake','ace_e','lac_e','fum_e','bio_yield']
    df = pd.read_excel(filename, sheet_name= 'Physiology', header=0, index_col=0).dropna(axis=0, how='all')
    df = pd.DataFrame(df.values, index=short_names, columns=df.columns).replace('-','NaN')
    for col in df.columns:
        df.loc[:,col] = df.loc[:,col].str.split(u' ± ')
    s = df.stack().apply(pd.Series)
    mean = s[0].unstack()
    sd = s[1].unstack()
    sd.columns = [x+'_sd' for x in sd.columns]
    df = pd.concat([mean,sd], axis=1)
    return df

class BiomassComposition(object):
    
    def __init__(self, organism='e_coli_fixed'):
        self.organism = organism
        
        if self.organism is not None:
            self.precursor_df = pd.read_excel(settings.DATA_DIR+'/meta_analysis.xls',
                                              sheet_name=self.organism + "_biomass",
                                              index_col=0, header=0).fillna(0)
            
            # a Series containing the growth dependence function for each sector
            self.growth_dependence = self.precursor_df['growth dependence'].apply(
                    lambda s: sympy.sympify(s))
            
            # a DataFrame with the required amount of each precursor in mmol
            # per 1 gram of that sector
            self.precursor_df = self.precursor_df.iloc[:, 1:]
    
    def GetSectorCorrections(self, growth_rate):
        """
            Returns a dataframe with the sector correction factors, which are
            the only growth-rate dependent part
        """
        mu = sympy.symbols('mu')
        
        if type(growth_rate) == float:
            growth_rate = ufloat(growth_rate, 0)
        
        sectors = {}
        for sect, func in self.growth_dependence.items():
            factor = func.evalf(subs={mu:growth_rate.nominal_value})
            uncert = np.abs(sympy.diff(func).evalf(subs={mu:growth_rate})) \
                * growth_rate.std_dev
            sectors[sect] = ufloat(factor, uncert)
        
        return pd.Series(sectors)
    
    def GetComposition(self, growth_rate, uncertainty=0.0):
        """
            Returns a Series containing values of each precursor in the biomass
            in units of mmol per gram of cell dry weight
        """
        if self.organism is None:
            return pd.Series()
        else:
            sector_correction_factors = self.GetSectorCorrections(growth_rate)
            biomass_df = self.precursor_df.transpose().multiply(sector_correction_factors)
            return biomass_df.transpose().apply(np.sum)

def adjust_biomass(model, new_biomass):
    """
        Adjust the biomass reaction in accordance with the growth rate
    """
    mets = []
    mapping = {'ATP':'atp_c = adp_c + pi_c',
               'G6P':'g6p_c',
               'PEP':'pep_c',
               'PYR':'pyr_c',
               'F6P':'g6p_c',
               'T3P':'dhap_c',
               'PGA':'3pg_c',
               'NADH':'nadh_c = nad_c',
               'NADPH':'nadph_c = nadp_c',
               'CO2':'co2_c',
               'R5P':'r5p_c',
               'E4P':'e4p_c',
               'OAA':'oaa_c',
               'AcCoA':'accoa_c = coa_c',
               'OGA':'akg_c'}
    
    for k,v in new_biomass.to_dict().items():
        if k == 'ATP': # 'ATP':'atp_c = adp_c + pi_c',
            mets += [('atp_c', v.nominal_value * -1)]
            mets += [('adp_c', v.nominal_value)]
            mets += [('pi_c', v.nominal_value)]
        elif k == 'NADH':
            mets += [('nadh_c', v.nominal_value * -1)]
            mets += [('nad_c', v.nominal_value)]
        elif k == 'NADPH':
            mets += [('nadph_c', v.nominal_value * -1)]
            mets += [('nadp_c', v.nominal_value)]
        elif k == 'AcCoA':
            mets += [('accoa_c', v.nominal_value * -1)]
            mets += [('coa_c', v.nominal_value)]    
        else:   
            mets += [(mapping[k], v.nominal_value * -1)]
        
    ### remove the old, replace by new
    model.remove_reactions(['Biomass_Ecoli_core_w_GAM'])
    Biomass_Ecoli_core_w_GAM = cobra.Reaction('Biomass_Ecoli_core_w_GAM') 
    d = {model.metabolites.get_by_id(met[0]): met[1] for met in mets}
    Biomass_Ecoli_core_w_GAM.add_metabolites(d)
    model.add_reactions([Biomass_Ecoli_core_w_GAM])
    model.objective = 'Biomass_Ecoli_core_w_GAM'
    
    return model


def main():
    """ Reconstruction of the model from Gerosa, 2015.
    Analysis flux imbalance with respect to the reported fluxes, and corrects
    for errors in the published supplementary data. 
    The reported flux data is mapped to the BiGG name space (as in iJO1366), 
    then assigned to individual reactions such that directionality matches the 
    extended core model to facilitate the projection of fluxes (using MOMA).
    DOI: http://dx.doi.org/10.1016/j.cels.2015.09.008
    """
    
    gerosa_balance_outfile = settings.CACHE_DIR+'/gerosa_imbalance.csv'
    gerosa_fluxes_outfile = settings.CACHE_DIR+'/gerosa_fluxes.csv'
    
    def Gerosa_model(gerosa_model_excel):
        """construct the gerosa model from the supplementary file
        """
        def add_mets(mets, reaction, sign): 
            for sub in mets.split(' + '):
                if ' ' in sub.strip():
                    c, m = sub.strip().split(' ')
                else:
                    c = 1
                    m = sub.strip()
                ## fixing some name space issues
                m = m.replace('-','__')
                if m == 'acon__C_c':
                    m = 'acon_C_c'
                m = cobra.Metabolite(m)
                reaction.add_metabolites({m:sign*float(c)})
            return reaction
            
        df = pd.read_excel(gerosa_model_excel, sheet_name= 'reactions', header=0, index_col=0)
    
        ## fixing name space by manual mapping, adopting iJO1366 name space
        gerosa_mapping = {'FUM_SEC':'FUMt2_2','SUCCt':'SUCCt2_2','PYRt2r':'PYRt2r',
                   'GlcnUpt':'GLCNt2r','GlycUpt':'GLYCt','G3PD5':'G3PD2','GLCpts':'GLCpts',
                   'FRUpts':'FRUpts','GALabc':'GALabc','ZWF':'G6PDH2r','ACONT1':'ACONTa','ACONT2':'ACONTb',
                   'SUCDH3':'SUCDi','MQO':'MOX','Growth_rate':'Biomass_Ecoli_core_w_GAM','LDH_D':'LDH_D',
                   'EX_glc_e':'EX_glc__D_e','EX_gln_L_e':'EX_gln__L_e','EX_glu_L_e':'EX_glu__L_e',
                   'EX_lac_D_e':'EX_lac__D_e','EX_mal_L_e':'EX_mal__L_e'}
    
        model = cobra.Model('gerosa_model')
        for rxn in df.index:
            name = rxn.replace('(','_').replace(')','')
            if name in gerosa_mapping:
                name = gerosa_mapping[name]
            reaction = cobra.Reaction(name)
            reaction.bounds = df.loc[rxn,['LB','UB']]
            formula = df.loc[rxn,'Formula'].replace('[','_').replace(']','')
            rev = ' <=> '; irrev = ' -> '
            if rev in formula:
                subs, prods = formula.split(rev)
            elif irrev in formula:
                subs, prods = formula.split(irrev)
            if subs:
                reaction = add_mets(subs, reaction, -1)
            if prods:
                reaction = add_mets(prods, reaction, 1)
            model.add_reactions([reaction])
    
        for rxn in model.reactions: # allow for secretion
            if rxn.boundary:
                model.reactions.get_by_id(rxn.id).upper_bound = 1000
                
        return model
    
    def extract_stoich(model):
        """extracts the stoichiometric matrix, outputs a dataframe"""
        columns = [rxn.id for rxn in model.reactions]
        index = [met.id for met in model.metabolites]
        S = pd.DataFrame(create_stoichiometric_matrix(model), index=index, columns=columns)
        S.columns.name = 'reaction'
        S.index.name = 'metabolite'
        return S
    
    def Flux_13C(filename):
        """Parsing and mapping of the fluxes to the proper name space"""
        rxn_mapping = {'GLCpts':'GLCpts , EX_glc__D_e',
                       'GlcnUpt':'GLCNt2r , EX_glcn_e',
                       'GALabc':'GALabc , EX_gal_e',
                       'PYRt2r':'PYRt2r , EX_pyr_e',
                       'GlycUpt':'GLYCt , EX_glyc_e',
                       'SUCCt':'SUCCt2_2 , EX_succ_e',
                       'ACS':'ACS , EX_ac_e',
                       'FRUpts':'FRUpts , EX_fru_e',
                       'FUM_SEC':'FUMt2_2 , EX_fum_e',
                       'LDH_D':'LDH_D , EX_lac__D_e',
                       'G3PD5':'G3PD2',
                       'ZWF':'G6PDH2r',
                       'ACONT':'ACONTa',
                       'ACONT2':'ACONTb',
                       'SUCDH3':'SUCDi',
                       'MQO':'MOX',
                       'Growth_rate':'Biomass_Ecoli_core_w_GAM',
                       }
    
        df = pd.read_excel(filename, sheet_name= 'Metabolic fluxes', header=0, index_col=0).dropna(axis=1, how='all')
        columns = list(df.values[:1][0,0:8]) + list([x+'_sd' for x in df.values[:1][0,0:8]])
        flux_rxns = list(df.iloc[0:-2,19].values) + [df.index[-1]]
        index = pd.Series([' '.join([rxn_mapping[rxn] if rxn in rxn_mapping else rxn for rxn in re.findall(r'\b\w+\b|\+|\-|\(|\)|,', rxns)]) for rxns in flux_rxns], index=df.index[1:])
        return pd.DataFrame(list(df.values[1:][0:,0:16]) , columns=columns, index=index)
    
    def Expand_flux_13C(flux_13C):
        """Expanding measured fluxes by assigning them to individual reactions.
        Additionally, several manual curations are made to account for missing 
        information and errors in the published supplementary data.
        """
        ### fix coa_c / accoa_c / ac_c balance imbalance
# =============================================================================
#         b_c = 3.7478
#         COA = - flux_13C.loc['( ACt2r , ACKr , PTAr ) + ACS , EX_ac_e',:] \
#             + flux_13C.loc['( ICL , MALS )',:] \
#             + flux_13C.loc['( CS , ACONTa , ACONTb )',:] - flux_13C.loc['( PDH )',:]\
#             + b_c * flux_13C.loc['Biomass_Ecoli_core_w_GAM',:] 
#         flux_13C.loc['( ACt2r , ACKr , PTAr ) + ACS , EX_ac_e',:] = flux_13C.loc['( ACt2r , ACKr , PTAr ) + ACS , EX_ac_e',:].add(COA[0:8])
# =============================================================================
    
        vector_index = [item for sublist in [re.findall(r'\b\w+\b', flux) for flux in flux_13C.index] for item in sublist]
        df = pd.DataFrame(index = vector_index, columns = flux_13C.columns)
        for flux in flux_13C.index:
            rxns = re.findall(r'\b\w+\b', flux)
            signs = [np.sign(x) > 0 for x in flux_13C.loc[flux,:]][:8] * 2 # mean sign determines assignment of sd to rxns
            
            ### futile cycles
            if flux == '( PYK - PPS )' or flux == '( PFK - FBP )':
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]  * signs
                df.loc[rxns[1],:] = -flux_13C.loc[flux,:] * [not i for i in signs]
            
            ### isoenzymes (parallel reactions)
            elif flux == '( MDH + MOX )' or flux == '( ME1 + ME2 )':
                df.loc[rxns[0],:] = flux_13C.loc[flux,:] * signs
                df.loc[rxns[1],:] = [0]*16
                
            ### rest from here is uptake systems
            elif flux == '( GLCpts , EX_glc__D_e )':
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]
                df.loc[rxns[1],:] = -flux_13C.loc[flux,:]
            elif flux == '( GLCNt2r , EX_glcn_e , GNK )':
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]
                df.loc[rxns[1],:] = -flux_13C.loc[flux,:]
                df.loc[rxns[2],:] = flux_13C.loc[flux,:]
            elif flux == '( GALabc , EX_gal_e , GALKr , UGLT , PGMT )': 
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]
                df.loc[rxns[1],:] = -flux_13C.loc[flux,:]
                df.loc[rxns[2],:] = flux_13C.loc[flux,:]
                df.loc[rxns[3],:] = flux_13C.loc[flux,:]
                df.loc[rxns[4],:] = flux_13C.loc[flux,:]
            elif flux == '( PYRt2r , EX_pyr_e )':
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]
                df.loc[rxns[1],:] = -flux_13C.loc[flux,:]
            ### we have determined that the glycerol uptake needs to be twice as 
            ### high as reported to match the measured fluxes.
            elif flux == '( GLYCt , EX_glyc_e , GLYK , G3PD2 )': 
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]  * 2
                df.loc[rxns[1],:] = -flux_13C.loc[flux,:] * 2
                df.loc[rxns[2],:] = flux_13C.loc[flux,:]  * 2
                df.loc[rxns[3],:] = flux_13C.loc[flux,:]  * 2
            elif flux == '( SUCCt2_2 , EX_succ_e )':
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]
                df.loc[rxns[1],:] = -flux_13C.loc[flux,:]
            elif flux == '( ACt2r , ACKr , PTAr ) + ACS , EX_ac_e':
                print(signs)
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]
                df.loc[rxns[1],:] = flux_13C.loc[flux,:]  * [not i for i in signs]
                df.loc[rxns[2],:] = -flux_13C.loc[flux,:] * [not i for i in signs]
                df.loc[rxns[3],:] = flux_13C.loc[flux,:]  * signs
                df.loc[rxns[4],:] = -flux_13C.loc[flux,:]
            elif flux == '( FRUpts , EX_fru_e , FRUK )': 
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]
                df.loc[rxns[1],:] = -flux_13C.loc[flux,:]
                df.loc[rxns[2],:] = flux_13C.loc[flux,:]
            elif flux == '( LDH_D , EX_lac__D_e , D_LACt2 )':
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]
                df.loc[rxns[1],:] = -flux_13C.loc[flux,:]
                df.loc[rxns[2],:] = flux_13C.loc[flux,:]
            elif flux == '( FUMt2_2 , EX_fum_e )': 
                df.loc[rxns[0],:] = flux_13C.loc[flux,:]
                df.loc[rxns[1],:] = -flux_13C.loc[flux,:]
            
            ### set the rest of the intracellular reactions
            else:
                for rxn in rxns: 
                    df.loc[rxn,:] = flux_13C.loc[flux,:]
        
        ### fix akg_c balance
# =============================================================================
#         b_c = 4.1182
#         GLUSy = df.loc['ICDHyr',:] - df.loc['AKGDH',:] + b_c * df.loc['Biomass_Ecoli_core_w_GAM',:]
#         df = df.append(pd.DataFrame(np.append(np.array(GLUSy[0:8]), (np.zeros(8))), index = df.columns, columns=['GLUSy']).T, )
#         df = df.append(pd.DataFrame(np.append(np.array(GLUSy[0:8]), (np.zeros(8))), index = df.columns, columns=['GLNS']).T, )
# =============================================================================
    
        ### flip signs for reactions that are reversed
        for r in df.index:
            if r in ['PGM','PGK','RPI','SUCOAS']:
                df.loc[r, :] = -df.loc[r, :]
    
        ### condition specific flipping for our extended core model, not gerosa
# =============================================================================
#         for cond in ['Acetate','Glycerol','Galactose']:
#             df.loc['ACKr',cond] = -df.loc['ACKr',cond]
#         
#         for cond in ['Pyruvate']:
#             df.loc['LDH_D',cond]       = -df.loc['LDH_D',cond]
#             df.loc['EX_lac__D_e',cond] = -df.loc['EX_lac__D_e',cond]
#             df.loc['D_LACt2',cond]     = -df.loc['D_LACt2',cond]
#         
#         for cond in ['Acetate']: # flipping: ( ACKr , PTAr ) + ACS
#             df.loc['ACKr',cond] = df.loc['ACS',cond]
#             df.loc['ACS',cond]  = -df.loc['PTAr',cond]
#             df.loc['PTAr',cond] = -df.loc['ACKr',cond]    
#     
#         for cond in ['Succinate']:
#             df.loc['EX_succ_e',cond] = df.loc['EX_succ_e',cond] + 2.25
#         
# =============================================================================
        return df
    
    def determine_imbalance(gerosa_model_excel, rxns_13C, physiology):
        """Analyse and correct imbalances in the Gerosa model"""
        gerosa_model = Gerosa_model(gerosa_model_excel)
        ## adjust biomass according to growth rate
        
        S = extract_stoich(gerosa_model)
        fluxes = rxns_13C.copy()
        df = pd.DataFrame(index = S.index, columns = settings.CONDITIONS)
        for condition in settings.CONDITIONS:
            ebc = BiomassComposition('e_coli')
            growth_rate = float(physiology.loc['growth_rate',condition])
            uncertainty = float(physiology.loc['growth_rate',condition+'_sd'])
            new_biomass = ebc.GetComposition(growth_rate, uncertainty)
            model = adjust_biomass(gerosa_model.copy(), new_biomass)        
            S = extract_stoich(model)
            
            fluxes_full = S.transpose().reset_index().join(fluxes[condition], on=S.columns.name)
            fluxes_full = fluxes_full.set_index('reaction')[condition].fillna(0.0)
            imbalance = pd.Series(data = np.dot(S, fluxes_full).round(1), index=S.index)
            COFACTORS = ['h2o_c', 'h_e', 'h_c', 'atp_c', 'adp_c', 'nadh_c', 'nad_c',
                         'nadph_c', 'nadp_c', 'q8_c', 'q8h2_c', 'co2_c', 'pi_c', 'amp_c',]
            imbalance[COFACTORS] = 0
            sol = imbalance[imbalance != 0]
            for met in sol.index:
                df.loc[met,condition] = sol[met]
        df = df.dropna(axis=0, how='all')#.fillna(0) 
        ## acetate imbalance on acetate due to condition-specific flipping acetate condition
        df.to_csv(gerosa_balance_outfile)
        

    gerosa_model_excel = settings.GEROSA_S4
    flux_13C = Flux_13C(settings.GEROSA_S2)
    rxns_13C = Expand_flux_13C(flux_13C)
    physiology = Physiology(settings.GEROSA_S2)    
    
    
    determine_imbalance(gerosa_model_excel, rxns_13C, physiology)
    print("Gerosa flux imbalance determined, can be found here: {}".format(gerosa_balance_outfile))
    
    rxns_13C.to_csv(gerosa_fluxes_outfile)
    print("Gerosa fluxes for flux projection can be found here: {}".format(gerosa_fluxes_outfile))
    
if __name__ == '__main__':
    main()
    
    
    
        
    