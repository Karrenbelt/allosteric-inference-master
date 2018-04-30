#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 16:37:19 2018

@author: Zarathustra
"""
import settings
import numpy as np
import pandas as pd
import sympy
import cobra
from uncertainties import ufloat

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

def adjust_biomass(model, growth_rate, uncertainty):
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
    
    
    ebc = BiomassComposition('e_coli')
    new_biomass = ebc.GetComposition(growth_rate, uncertainty)    
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

            