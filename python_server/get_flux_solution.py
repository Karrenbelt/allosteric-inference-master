#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 12:34:09 2018

@author: Zarathustra
"""

import settings
import pandas as pd
import numpy as np
import sympy
import cobra
import cobra.util.solver as sutil
from optlang.symbolics import Zero
import matplotlib.pyplot as plt
from uncertainties import ufloat

# =============================================================================
# def nullify_objectives(model):
#     """
#         Nullifies all the objectives in the constraint-based model
#     """
#     [setattr(x, 'objective_coefficient', 0) for x in model.reactions if x.objective_coefficient != 0]
# 
# def get_medium(model):
#     """
#         returns the medium composition in the constraint-based model
#     """
#     return pd.DataFrame([model.medium[met] for met in model.medium], index=[met for met in model.medium], columns=['mM'])
# 
# def constrain_ATPM(model):
#     """
#         maximixes ATP production and sets it as a constaint
#     """
#     current_objectives = [(rxn.id, rxn.objective_coefficient)  for rxn in model.reactions if model.reactions.get_by_id(rxn.id).objective_coefficient!=0]
#     nullify_objectives(model)
#     model.objective = "ATPM"
#     model.reactions.get_by_id('ATPM').upper_bound = model.optimize().objective_value
#     nullify_objectives(model)
#     for rxn in current_objectives:
#         model.reactions.get_by_id(rxn[0]).objective_coefficient = rxn[1]
# 
# def FBA(model, condition, rxns_13C, method='optgp', n_samples=1):
#     """
#         Obtain a FBA solution using the reported uptake rate only
#     """
#     uptake = settings.COND2BIGG[condition]
#     model.reactions.get_by_id(uptake).lower_bound = rxns_13C.loc[uptake, condition]
#     solution = model.optimize().fluxes
#     solution.to_json(settings.CACHE_DIR+'/FBA_'+condition+'.json')
# #     solution = cobra.flux_analysis.sample(model, n_samples, method=method) # method = 'arch' 
#     return solution.to_frame()
# 
# def FVA(model, condition, rxns_13C):
#     """
#         Flux variance analysis returns min, max flux range given 99% optimality
#     """
#     uptake = settings.COND2BIGG[condition]
#     model.reactions.get_by_id(uptake).lower_bound = rxns_13C.loc[uptake,condition]
#     solution = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=.99)
# #     solution.to_json(output+'/FVA_'+condition+'.json') # not suitable for Escher, need to sample the space first?
#     return solution # already dataframed    
# 
# def pFBA(model, condition, rxns_13C):
#     """
#         Parsimonious FBA returns optimal solution with the minimal sum of fluxes
#     """
#     uptake = settings.COND2BIGG[condition]
#     model.reactions.get_by_id(uptake).lower_bound = rxns_13C.loc[uptake,condition]
#     solution = cobra.flux_analysis.parsimonious.pfba(model, objective="Biomass_Ecoli_core_w_GAM").fluxes.to_frame()
#     solution.to_json(settings.CACHE_DIR+'/pFBA_'+condition+'.json') 
#     return solution
# =============================================================================

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

def project_fluxes(model, condition, rxns_13C):
    """
        Using MOMA to find a flux solution as close to the published data as possible
    """
    #print(settings.CONDITIONS.keys())
    exchange_2_ID = settings.CONDITIONS
    ### MOMA: obtain flux solution that minimally deviates from 13C flux measurements 
    prob = model.problem
    v = prob.Variable("moma_old_objective")
    c = prob.Constraint(model.solver.objective.expression - v, lb=0.0, ub=0.0, name="moma_old_objective_constraint")
    to_add = [v, c]
    new_obj = Zero
    
    linear = True
    fluxes_df = pd.DataFrame(index=rxns_13C.index, columns=['measured','sigma','predicted','lb','ub'])
    gen = [r for r in model.reactions if r.id in rxns_13C.index]
    
    for r in gen:    
        flux = rxns_13C.loc[r.id,condition]
        
        if r.id in [exchange_2_ID[condition]]:
            model.reactions.get_by_id(exchange_2_ID[condition]).lower_bound = flux
        elif r.id.startswith('EX_'):
            ### in the succinate condition there is also fumerate uptake:
            ### must have been secreted by the cells first
            if condition == 'Succinate' and r.id == 'EX_fum_e':
                model.reactions.get_by_id('EX_fum_e').lower_bound = flux 
            else:
                flux = abs(flux) ## secretion
        
        fluxes_df.loc[r.id,'measured'] = flux
        fluxes_df.loc[r.id,'sigma'] = rxns_13C.loc[r.id,condition+'_sd']
        
        if linear:
            model.solver = 'glpk'
            components = sutil.add_absolute_expression(model, r.flux_expression, name="moma_dist_" + r.id, difference=flux, add=False)
            to_add.extend(components)
            new_obj += components.variable
        else:
            model.solver = 'cplex' # not working properly
            dist = prob.Variable("moma_dist_" + r.id)
            const = prob.Constraint(r.flux_expression - dist, lb=flux, ub=flux, name="moma_constraint_" + r.id)
            to_add.extend([dist, const])
            new_obj += dist**2
            
    model.add_cons_vars(to_add)
    model.objective = prob.Objective(new_obj, direction='min')

    solution = model.optimize().fluxes
    solution.to_json(settings.CACHE_DIR+'/FBA_'+condition+'.json')
    
    ### make scatter plot of predicted VS measured
    for r in gen:
        fluxes_df.loc[r.id,'predicted'] = solution.loc[r.id]
        fluxes_df.loc[r.id,'lb'] = model.reactions.get_by_id(r.id).lower_bound
        fluxes_df.loc[r.id,'ub'] = model.reactions.get_by_id(r.id).upper_bound
    fig, axs = plt.subplots(1, 2, figsize=(14,6))
    axs[0].plot([fluxes_df['measured'].min(), fluxes_df['measured'].max()], [fluxes_df['predicted'].min(), fluxes_df['predicted'].max()], 'k--', alpha=0.3, linewidth=0.5)
    plot = fluxes_df.plot(kind='scatter', x=['measured'], y=['predicted'], xerr=fluxes_df.loc[:,'sigma'], 
          title = condition, ax=axs[0], linewidth=0, s=10, color=(0.7,0.2,0.5))
    
    ### annotate reactions
    for r in gen:
        xy = fluxes_df.loc[r.id, ['measured', 'predicted']]
        sign = np.sign(xy['measured']) == np.sign(xy['predicted'])
        if not sign:
            axs[0].annotate(r.id, xy, xytext=(10,-5), textcoords='offset points',
                        family='sans-serif', fontsize=10, color='darkslategrey')
    
    ### histogram of residuals
    fluxes_df['diff'] = fluxes_df[['measured']].sub(fluxes_df['predicted'], axis=0)
    residual = fluxes_df.loc[:, 'diff']
    abs(residual[residual!=0]).plot(kind='bar', ax=axs[1], color=(0.7,0.2,0.5))
    axs[1].set_xlabel('residual [mmol/gCDW/h]')        
    fig.savefig(settings.CACHE_DIR+'/MOMA_'+condition+'.pdf')
    
    return solution.to_frame()#, fluxes_df


def condition_specific_flip(df):
    
# =============================================================================
#     for cond in ['Acetate','Glycerol','Galactose']:
#         df.loc['ACKr',cond] = -df.loc['ACKr',cond]
#     
#     for cond in ['Pyruvate']:
#         df.loc['LDH_D',cond]       = -df.loc['LDH_D',cond]
#         df.loc['EX_lac__D_e',cond] = -df.loc['EX_lac__D_e',cond]
#         df.loc['D_LACt2',cond]     = -df.loc['D_LACt2',cond]
#     
#     for cond in ['Acetate']:  # flipping: ( ACKr , PTAr ) + ACS
#         df.loc['ACKr',cond] = df.loc['ACS',cond]
#         df.loc['ACS',cond]  = -df.loc['PTAr',cond]
#         df.loc['PTAr',cond] = -df.loc['ACKr',cond]    
#          
#     for cond in ['Succinate']:
#         df.loc['EX_succ_e',cond] = df.loc['EX_succ_e',cond] + 2.25
# =============================================================================

    return df

def main():
    """
        Obtain a flux solution for the extended core model, for the different
        conditions (carbon sources) studied in Gerosa (2015)
    """
    rxns_13C = pd.read_csv(settings.ECOLI_GEROSA_FLUX, index_col=0)
    rxns_13C = condition_specific_flip(rxns_13C)
    
    extended_core = cobra.io.read_sbml_model(settings.ECOLI_EXCORE_FNAME)
    physiology = Physiology(settings.GEROSA_S2)
    projected_flux_outfile = settings.CACHE_DIR+'/flux_solutions.p'
    
    ## project fluxes on the gerosa model
    flux_solution = dict()
    for condition in settings.CONDITIONS.keys():
        print(condition)
        if True: # mode == 1 (flux projection)
            ebc = BiomassComposition('e_coli')
            growth_rate = float(physiology.loc['growth_rate',condition])
            uncertainty = float(physiology.loc['growth_rate',condition+'_sd'])
            new_biomass = ebc.GetComposition(growth_rate, uncertainty)
            model = adjust_biomass(extended_core.copy(), new_biomass)
            flux_solution[condition] = project_fluxes(model, condition, rxns_13C)
        else:
            flux_solution[condition] = FBA(extended_core.copy(), condition, rxns_13C)
    
    pd.Panel(flux_solution).to_pickle(projected_flux_outfile)
    print('A flux solution generated, writen here: {}'.format(projected_flux_outfile))

if __name__ == '__main__':
    main()