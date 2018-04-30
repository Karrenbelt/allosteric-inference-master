#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 12:00:30 2018

@author: Zarathustra
"""

import settings, cobra, re
import pandas as pd
import numpy as np
from cobra.util import create_stoichiometric_matrix
import cobra.util.solver as sutil
from optlang.symbolics import Zero
import matplotlib.pyplot as plt
from biomass import adjust_biomass

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
    df.to_csv(settings.CACHE_DIR+'/physiology.csv')
    print('Gerosa physiological data writen here: {}'.format(settings.CACHE_DIR+'/gerosa_physiology.csv'))
    return df

def Gerosa_model(gerosa_model_excel):
    """construct the gerosa model from the supplementary file
       also adopts the name space of iJO1366 (genome-scale model)
    """
    def add_mets(mets, reaction, sign): 
        for sub in mets.split(' + '):
            if ' ' in sub.strip():
                c, m = sub.strip().split(' ')
            else:
                c = 1
                m = sub.strip()

            m = m.replace('-','__')
            if m == 'acon__C_c':
                m = 'acon_C_c'
            m = cobra.Metabolite(m)
            reaction.add_metabolites({m:sign*float(c)})
        return reaction
        
    df = pd.read_excel(gerosa_model_excel, sheet_name='reactions', header=0, index_col=0)

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

def flux_13C(filename):
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

def expand_flux_13C(flux_13C):
    """Expanding measured fluxes by assigning them to individual reactions.
    Additionally, several manual curations are made to account for missing 
    information and errors in the published supplementary data.
    """
    ### fix coa_c / accoa_c / ac_c balance imbalance
    b_c = 3.7478
    COA = - flux_13C.loc['( ACt2r , ACKpr , PTAr ) + ACS , EX_ac_e',:] \
        + flux_13C.loc['( ICL , MALS )',:] \
        + flux_13C.loc['( CS , ACONTa , ACONTb )',:] - flux_13C.loc['( PDH )',:]\
        + b_c * flux_13C.loc['Biomass_Ecoli_core_w_GAM',:] 
    flux_13C.loc['( ACt2r , ACKr , PTAr ) + ACS , EX_ac_e',:] = flux_13C.loc['( ACt2r , ACKr , PTAr ) + ACS , EX_ac_e',:].add(COA[0:8])

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

    return df

def determine_imbalance(gerosa_model_excel, rxns_13C, physiology, outfile):
    """Analyse and correct imbalances in the Gerosa model"""
    gerosa_model = Gerosa_model(gerosa_model_excel)
    ## adjust biomass according to growth rate
    
    S = extract_stoich(gerosa_model)
    fluxes = rxns_13C.copy()
    df = pd.DataFrame(index = S.index, columns = settings.CONDITIONS)
    for condition in settings.CONDITIONS:
        growth_rate = float(physiology.loc['growth_rate',condition])
        uncertainty = float(physiology.loc['growth_rate',condition+'_sd'])
        model = adjust_biomass(gerosa_model.copy(), growth_rate, uncertainty)        
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
    df.to_csv(outfile)
    
    
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
    
    #
# =============================================================================
#     S = extract_stoich(model)
#     imbalance = pd.Series(data = np.dot(S, solution), index=S.index)
#     print(imbalance[abs(imbalance) > 1e-10])
#     solution.to_json(settings.CACHE_DIR+'/FBA_'+condition+'.json')
# =============================================================================
    
    # make scatter plot of predicted VS measured
    for r in gen:
        fluxes_df.loc[r.id,'predicted'] = solution.loc[r.id]
        fluxes_df.loc[r.id,'lb'] = model.reactions.get_by_id(r.id).lower_bound
        fluxes_df.loc[r.id,'ub'] = model.reactions.get_by_id(r.id).upper_bound
    fig, axs = plt.subplots(1, 2, figsize=(14,6))
    axs[0].plot([fluxes_df['measured'].min(), fluxes_df['measured'].max()], [fluxes_df['predicted'].min(), fluxes_df['predicted'].max()], 'k--', alpha=0.3, linewidth=0.5)
    fluxes_df.plot(kind='scatter', x=['measured'], y=['predicted'], xerr=fluxes_df.loc[:,'sigma'], 
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
    
    for cond in ['Acetate','Glycerol','Galactose']:
        df.loc['ACKr',cond] = -df.loc['ACKr',cond]
    
    for cond in ['Pyruvate']:
        df.loc['LDH_D',cond]       = -df.loc['LDH_D',cond]
        df.loc['EX_lac__D_e',cond] = -df.loc['EX_lac__D_e',cond]
        df.loc['D_LACt2',cond]     = -df.loc['D_LACt2',cond]
    
    for cond in ['Acetate']:  # flipping: ( ACKr , PTAr ) + ACS
        df.loc['ACKr',cond] = df.loc['ACS',cond]
        df.loc['ACS',cond]  = -df.loc['PTAr',cond]
        df.loc['PTAr',cond] = -df.loc['ACKr',cond]    
         
    for cond in ['Succinate']:
        df.loc['EX_succ_e',cond] = df.loc['EX_succ_e',cond] + 2.25

    return df

def main():
    """
        Obtain a flux solution for the extended core model, for the different
        conditions (carbon sources) studied in Gerosa (2015)
    """

    ## inputs
    gerosa_model_excel = settings.GEROSA_S4
    extended_core = cobra.io.read_sbml_model(settings.ECOLI_EXCORE_FNAME)
    fluxes = flux_13C(settings.GEROSA_S2)
    flux_per_rxn = expand_flux_13C(fluxes)
    physiology = Physiology(settings.GEROSA_S2)    
    
    ## output
    gerosa_balance_outfile = settings.CACHE_DIR+'/gerosa_imbalance.csv'
    gerosa_fluxes_outfile = settings.CACHE_DIR+'/gerosa_fluxes.csv'
    projected_flux_outfile = settings.CACHE_DIR+'/flux_solutions.p'   

    
    determine_imbalance(gerosa_model_excel, flux_per_rxn, physiology, gerosa_balance_outfile)
    print("Gerosa flux imbalance determined, can be found here: {}".format(gerosa_balance_outfile))
    
    flux_per_rxn.to_csv(gerosa_fluxes_outfile)
    print("Gerosa fluxes for flux projection can be found here: {}".format(gerosa_fluxes_outfile))
        
    
    flux_per_rxn = condition_specific_flip(flux_per_rxn)

    ## project fluxes on the gerosa model
    flux_solution = dict()
    for condition in settings.CONDITIONS.keys():
        if True: # mode == 1 (flux projection)
# =============================================================================
#             growth_rate = float(physiology.loc['growth_rate',condition])
#             uncertainty = float(physiology.loc['growth_rate',condition+'_sd'])
#             extended_core = adjust_biomass(extended_core.copy(), growth_rate, uncertainty)
# =============================================================================
            flux_solution[condition] = project_fluxes(extended_core.copy(), condition, flux_per_rxn)
        else:
            flux_solution[condition] = FBA(extended_core.copy(), condition, flux_per_rxn)
    
    pd.Panel(flux_solution).to_pickle(projected_flux_outfile)
    print('A flux solution generated, writen here: {}'.format(projected_flux_outfile))

if __name__ == '__main__':
    main()    
    
    
    
        
    