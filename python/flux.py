import settings
import os, re, cobra
import pandas as pd
import numpy as np

# def nullify_objectives(model):
#     [setattr(x, 'objective_coefficient', 0) for x in model.reactions if x.objective_coefficient != 0]

# def get_truncnorm(mean=0, sd=1, lb=0, ub=10):
#     return truncnorm((lb - mean) / sd, (ub - mean) / sd, loc=mean, scale=sd)

# def get_medium(model):
#     return pd.DataFrame([model.medium[met] for met in model.medium], index=[met for met in model.medium], columns=['mM'])

# def constrain_ATPM(model): # currently not used
#     current_objectives = [(rxn.id, rxn.objective_coefficient)  for rxn in model.reactions if model.reactions.get_by_id(rxn.id).objective_coefficient!=0]
#     nullify_objectives(model)
#     model.objective = "ATPM"
#     model.reactions.get_by_id('ATPM').upper_bound = model.optimize().objective_value
#     nullify_objectives(model)
#     for rxn in current_objectives:
#         model.reactions.get_by_id(rxn[0]).objective_coefficient = rxn[1]

def FBA(model, condition, rxns_13C, method='optgp', n_samples=1):
    uptake = settings.COND2BIGG[condition]
    model.reactions.get_by_id(uptake).lower_bound = rxns_13C.loc[uptake, condition]
    solution = model.optimize().fluxes
    solution.to_json(settings.CACHE_DIR+'/FBA_'+condition+'.json')
#     solution = cobra.flux_analysis.sample(model, n_samples, method=method) # method = 'arch' 
    return solution.to_frame()

# def FVA(model, condition, rxns_13C):
#     uptake = settings.COND2BIGG[condition]
#     model.reactions.get_by_id(uptake).lower_bound = rxns_13C.loc[uptake,condition]
#     solution = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=.99)
# #     solution.to_json(output+'/FVA_'+condition+'.json') # not suitable for Escher, need to sample the space first?
    # return solution # already dataframed    

# def pFBA(model, condition, rxns_13C):
#     uptake = settings.COND2BIGG[condition]
#     model.reactions.get_by_id(uptake).lower_bound = rxns_13C.loc[uptake,condition]
#     solution = cobra.flux_analysis.parsimonious.pfba(model, objective="Biomass_Ecoli_core_w_GAM").fluxes.to_frame()
#     solution.to_json(settings.CACHE_DIR+'/pFBA_'+condition+'.json') # not suitable for Escher, need to sample the space first?
#     return solution


rxns_13C = pd.read_csv(settings.CACHE_DIR+'/gerosa_fluxes_per_rxn.csv', index_col=0)
extended_core = cobra.io.read_sbml_model(settings.ECOLI_EXCORE_FNAME)

flux_solution = dict()
for condition in settings.CONDITIONS:
    flux_solution[condition] = FBA(extended_core.copy(), condition, rxns_13C)

pd.Panel(flux_solution).to_pickle(settings.CACHE_DIR+'/flux_solutions.p')



