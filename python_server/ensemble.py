#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 16:42:44 2018

@author: Zarathustra
"""

import tellurium as te # needs to be imported before cobra, otherwise: "Segmentation fault: 11"
import cobra
import settings, sys, os, json, itertools, pickle
from cobra.util import create_stoichiometric_matrix
import pandas as pd
import numpy as np
from scipy.stats import truncnorm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def remove_mets_rxns(S, rxns=[], mets=[]):
    """
        Removal of reactions and metabolites from the stoichiometric matrix
    """
    ## remove metabolites
    if not set(mets).issubset(set(S.index)):
        no_match = list(set(mets).difference(set(S.index)))
        raise ValueError('Metabolites %s do not occur in the Stoichiometrix matrix' % no_match)
    else:
        S_sub = S.drop(mets)
        if not S_sub.columns.equals(S_sub.dropna(axis=1, how='all').columns):
            raise AssertionError('removal of metabolites results in empty reactions')

    ## remove reactions        
    if not set(rxns).issubset(set(S.columns)):
        no_match = list(set(rxns).difference(set(S.columns)))
        raise ValueError('Reactions %s do not occur in the Stoichiometrix matrix' % no_match)
    else:
        S_sub = S_sub.drop(rxns, axis=1)
        if not (S_sub.T != 0).any().all():
            raise AssertionError('removal of reactions results in removal of metabolites')     
    return S_sub

def Delta_gibbs(filename):
    mapping = {'ZWF':'G6PDH2r','ACONT1':'ACONTa','ACONT2':'ACONTb','SUCDH3':'SUCDi','MQO':'MOX'}
    df = pd.read_excel(filename, sheet_name= 5, header=2, index_col=0).dropna(axis=1, how='all')
    columns = [str(x) for x in df.columns[:8].values] + [str(x)+'_sd' for x in df.columns[:8].values]
    index = [mapping[x] if x in mapping else x for x in df.index[:59]]
    return pd.DataFrame(df.values[:59,0:16], columns = columns, index = index).dropna(axis=0, how='all')

class MakeModel:
    """
        Generates a kinetic model from a stoichiometric matrix
    """
    def __init__(self, S):
    
        ## change name space metabolites and reactions
        new_index = ['S'+str(i) for i in range(len(S.index))]
        new_columns = ['R'+str(i) for i in range(len(S.columns))]
        self.mapping = dict(zip(list(S.index)+list(S.columns), new_index+new_columns))
        self.mapping.update({v: k for k, v in self.mapping.items()})
        S.columns = new_columns
        S.index = new_index
        
        self.S        = S
        self.P        = dict().fromkeys([x+'_0' for x in S.index])
        self.X0       = dict.fromkeys(S.index)
        self.rxns     = dict.fromkeys(S.columns)
        self.rates    = dict.fromkeys(S.columns)
        self.vmax     = dict()
        self.exchange = [x for x in S.index if self.mapping[x].endswith('_e')]
        self.tag      = None
        self.antimony = None
        
        ## add the reaction formula
        for rxn in self.rxns:
            subs, prods = [],[]
            mets = [str(self.S[rxn][i])+' '+i for i in (self.S[rxn]!=0)[self.S[rxn]!=0].index]
            for met in mets:
                c, m = met.split(' ')        
                coeff = float(c)
                new_m = str(abs(coeff))+' '+ m if (abs(coeff) != 1) else m
                if coeff < 0: 
                    subs += [new_m]
                else:
                    prods += [new_m]
            self.rxns[rxn] = ' + '.join(subs)+' -> '+' + '.join(prods)

    def add_rate(self, rxn, law='rev_MM'):
        f, b, MAp, MAs = ([] for i in range(4))
        mets = [str(self.S[rxn][i])+' '+i for i in (self.S[rxn]!=0)[self.S[rxn]!=0].index]
        E = 'E_'+rxn; kcat = 'kcat_'+rxn; Keq = 'Keq_'+rxn; v = 'v_'+rxn
        Vmaxf = E+'_'+kcat
        Vmaxb = 'Vmaxb_'+rxn
        
        self.P[v] = None; # self.P[Vmaxf] = None; # for steady-state parameterization
#         self.P[kcat] = None; self.P[E] = None; # for parameter balancing sbml file output
    
        for met in mets:
            c, m = met.split(' ')
            Km = '_'.join(['Km',rxn,m])
            coeff = float(c)
            new_m = m + ' ^ ' + str(abs(coeff)) if (abs(coeff) != 1) else m
            if coeff < 0:
                if law in ['irrev_MM', 'rev_MM', 'Haldane']:
                    self.P[Km] = None
                    f += [new_m + ' / ' + Km]
                MAs += [new_m]
            else:
                if law in ['rev_MM', 'Haldane']:
                    self.P[Km] = None
                    b += [new_m + ' / ' + Km]
                MAp += [new_m]
        
        # a boundary reaction (with a constant flux) # if using S instead of S_sub
#         if self.mapping[rxn].startswith('EX_'):
#             rate = Vmaxf
#             vmax = v
            
        if self.mapping[rxn] == settings.BIOMASS_ID:
            f = []
            for met in mets:
                c, m = met.split(' ')
                Km = '_'.join(['Km',rxn,m])
                self.P[Km] = None
                coeff = float(c)
                new_m = m
                if coeff < 0:
                    f += [Km + ' / ' + new_m]
                rate = '%s / ( 1 + %s )' % (Vmaxf, ' + '.join(f))
                vmax = '%s * ( 1 + %s )' % (v, ' + '.join(f))
#                rate = '%s * %s / ( 1 + %s )' % (kcat, E, ' + '.join(f))
        
        # irreversible mass action
        elif law == 'irrev_MA': 
            rate = '%s * %s' % (Vmaxf, ' * '.join(MAs))
            vmax = '%s / %s' % (v, ' * '.join(MAs))

        # irreversible Michaelis-Menten
        elif law == 'irrev_MM':
            num  = ' * '.join(f)
            den  = '( 1 + %s )'       % (' + '.join(f))
            rate = '( %s * %s / %s )' % (Vmaxf, num, den)
            vmax = '( %s * %s / %s )' % (v, den, num)

        # reversible mass action
        elif law == 'rev_MA': 
            rev = 'rev_'+rxn
            self.P[rev] = None
            v_f   = '%s * ( 1 + %s )'   % (v, rev)
            v_b   = '%s * %s'           % (v, rev)
            vmax  = '%s / ( %s )'       % (v_f, ' * '.join(MAs))
            vmaxb = '%s / ( %s )'       % (v_b, ' * '.join(MAp))
            rate  = '%s * %s - %s * %s' % (Vmaxf, ' * '.join(MAs), Vmaxb, ' * '.join(MAp))
            self.vmax[Vmaxb] = vmaxb
#             self.P[Vmaxb] = None
        
        # reversible Michaelis-Menten
        elif law == 'rev_MM': 
            rev = 'rev_'+rxn
            self.P[rev] = None
            num  = '( %s * %s - %s * %s )'  % (Vmaxf, ' * '.join(f), Vmaxb, ' * '.join(b))
            den  = '( 1 + %s + %s )'        % (' + '.join(f), ' + '.join(b))
            v_f   = '%s * ( 1 + %s )'       % (v, rev)
            v_b   = '%s * %s'               % (v, rev)
            vmax  = '%s * %s / ( %s )'      % (v_f, den, ' * '.join(f))
            vmaxb = '%s * %s / ( %s )'      % (v_b, den, ' * '.join(b))
            rate  = '%s / %s'               % (num, den)
            self.vmax[Vmaxb] = vmaxb
#             self.P[Vmaxb] = None
        
        # reversible Michaelis-Menten with haldane substitution    
        elif law == 'Haldane': 
            self.P[Keq] = None;
            rev  = '( 1 - ( %s / %s ) / %s )'  % (' * '.join(MAp), ' * '.join(MAs), Keq)
            num  = ' * '.join(f)
            den  = '( 1 + %s + %s )'           % (' + '.join(f), ' + '.join(b))
            rate = '( %s * %s * %s / %s )'     % (Vmaxf, rev, num, den)
            vmax = '( %s / ( %s * %s / %s ) )' % (v, rev, num, den)
#             rate = '( %s * %s * %s * %s / %s )'     % (kcat, E, rev, num, den)
        
        self.rates[rxn] = rate
        self.vmax[Vmaxf] = vmax

    def add_regulation(self, regulators):
        self.tag = regulators
        for i in range(0,len(regulators[0])):
            rxn = self.mapping[regulators[0][i]]
            met = self.mapping[regulators[1][i]]
            alpha = 'a_'+rxn+'_'+met
            self.P[alpha] = None
            self.rates[rxn] += ' * ( %s / %s_0 ) ^ %s' % (met, met, alpha)
    
    def parameterize(self, flux_solution): 
        
#         meanLink = pd.read_excel(settings.CACHE_DIR+'/mean_Link.xlsx', index_col=0) 
#         meanLink = pd.DataFrame(meanLink)
#         sdLink = pd.read_excel(settings.CACHE_DIR+'/sd_Link.xlsx', index_col=0) 
#         sdLink = pd.DataFrame(sdLink)  ## some are nan
#         K = pd.read_csv(settings.CACHE_DIR+'/parameters', sep= '\t', index_col=0) 
#         K = pd.DataFrame(K)

        def get_truncnorm(mean=0, sd=1, lb=0, ub=10):
            return truncnorm((lb - mean) / sd, (ub - mean) / sd, loc=mean, scale=sd)
        
        for k,v in self.P.items():
            if k in self.vmax.keys():
                pass
            elif self.mapping.get(k[2:]) in flux_solution.index: ## 'v_'+rxn
                flux_value = float(flux_solution.loc[self.mapping[k[2:]]])
                self.P[k] = str(settings.CONVERT_UNITS(flux_value))
#             elif k.lower() in K.index: ## assigning known Km- and Keq-values (case insensitive)
#                 l = k.lower()
#                 if math.isnan(K.at[l, '!Std']):
#                     self.P[k] = str(get_truncnorm(mean=K.at[l, '!Mean'], sd=1, lb=0.1, ub=10).rvs())
#                 else:
#                     self.P[k] = str(get_truncnorm(mean=K.at[l, '!Mean'], sd=K.at[l, '!Std'], lb=0.1, ub=10).rvs())
            elif k.startswith('rev_'): # sample reversibility
                self.P[k] = str(np.random.uniform(0.2, 0.8))
            else:
                self.P[k] = str(get_truncnorm(mean=1, sd=1, lb=0.1, ub=10).rvs())
                
        for k,v in self.X0.items(): 
#             if k in meanLink.index:
#                 if math.isnan(sdLink.at[k, '10']):
#                     val = str(get_truncnorm(mean=meanLink.at[k, '10'], sd=1, lb=0, ub=10).rvs())
#                 else:
#                     val = str(get_truncnorm(mean=meanLink.at[k, '10'], sd=sdLink.at[k, '10'], lb=0, ub=10).rvs())
#             else:
            val = str(get_truncnorm(mean=1, sd=1, lb=0.1, ub=3).rvs())
            self.X0[k] = val
            self.P[k+'_0'] = val
            
        for k,v in self.vmax.items():
            self.vmax[k] = ' '.join([x if x in self.X0 else x for x in v.split(' ')]) 
                
    def to_antimony(self):
    
        def convert_string(string):
            return ' '.join([self.mapping[x] if x in self.mapping else x for x in string.split(' ')])
        
        antimony_str = ''''''
        antimony_str += '// model *coli_core()'

        antimony_str += '\n\n// Reactions:'       
        for rxn, formula in self.rxns.items():
            formula = ' '.join(['$ '+x if x in self.exchange else x for x in formula.split(' ')])
            antimony_str += '\n'+rxn+': '+formula+'; '+self.rates[rxn] + ';'
        
        antimony_str += '\n\n// Initial conditions:\n'
        for met, value in self.X0.items():
            antimony_str += ''+met+' = '+value+'; '

        antimony_str += '\n\n// Parameters:\n'
        for par, value in self.P.items():
            if par not in self.vmax.keys():
                antimony_str += ''+par+' = '+value+'; '

        antimony_str += '\n\n// Derived parameters:'
        for rxn, vmax in self.vmax.items():
            vmax = calculate_vmax(self, vmax)
            self.P[rxn] = vmax
            antimony_str += '\n '+rxn+' = '+vmax+'; '
            
        self.antimony = antimony_str
        
    def write_SBML_mapping(self):
        r = te.loada(self.antimony)
        r.exportToSBML(settings.CACHE_DIR+'/extended_core_kinetic.xml', current=False)
        json.dump(self.mapping, open(settings.CACHE_DIR+"/mapping.json",'w'))
        species_and_parameters = {**self.P, **self.X0}
        json.dump(species_and_parameters, open(settings.CACHE_DIR+"/model_species_and_parameters.json",'w'))
        #w = csv.writer(open(settings.CACHE_DIR+"/mapping.csv", "w"))
        #for key, val in self.mapping.items(): 
        #    w.writerow([key, val])
            
def calculate_vmax(model, vmax):
    mapped_equation = [model.X0[x] if x in model.X0 else model.P[x] if x in model.P else x for x in vmax.split(' ')]
    return str(eval(' '.join(mapped_equation).replace('^','**')))

def perturb(model, condition, new_condition): # to glucose 
    """
        Adds an event to the antimony model, simulating a carbon source switch
    """
    R1 = model.mapping[settings.UPTAKE[condition]]
    R2 = model.mapping[settings.UPTAKE[new_condition]]
    VmaxR1 = 'E_'+R1+'_kcat_'+R1
    VmaxR2 = 'E_'+R2+'_kcat_'+R2
    uptake_rate = pd.read_csv(settings.ECOLI_GEROSA_PHYS, index_col=0).loc['uptake',new_condition]
    model.P['v_'+R2] = str(settings.CONVERT_UNITS(uptake_rate))
    
    value = calculate_vmax(model, model.vmax[VmaxR2])
    model.antimony += '\n\n //Events: '
    model.antimony += '\nat (time > 10): %s = %s;' % (VmaxR1, '0')  
    model.antimony += '\nat (time > 10): %s = %s;' % (VmaxR2, value)
    model.antimony += '\nat (time > 40): %s = %s;' % (VmaxR1, model.P[VmaxR1]) 
    model.antimony += '\nat (time > 40): %s = %s;' % (VmaxR2, model.P[VmaxR2])
    
    return model

def allostery(M,R,N):
    """
        generates all N-wise allosteric interaction topologies
        args:
            M (list) - metabolite names
            R (list) - reaction names
            N (int)  - N-wise combinations
        returns:
            a list of tuples of tuples [(rxn_tuple), (met_tuple)]
    """
    rxns = itertools.combinations(R, N)
    mets = itertools.product(M, repeat=N)
    allo = list(itertools.product(rxns, mets))
    return allo

def simulate_antimony_model(model, plot=False, save_results=True, time_points=None):
    """
        Simulate the model
        args:
            model        - model object (from MakeModel, parameterized)
            plot         - True or False (default: False)
            save_results - True or False (default: True)
            time_points  - Save only specified time points (default: None)
        return:
            simulated time series in a pandas DataFrame
    """
    r = te.loada(model.antimony)
    r.timeCourseSelections= ['time'] + r.getFloatingSpeciesIds() + r.getBoundarySpeciesIds()
    
    try:
        sim_data = r.simulate(0,70,100)
    except RuntimeError:
        return pd.DataFrame() 
    
    if plot:
        plt.figure(figsize=(15,10))
        new = [model.mapping[x] for x in r.getFloatingSpeciesIds()+ r.getBoundarySpeciesIds()]
        te.plotArray(sim_data, ylabel='Concentration (mM)', xlabel='Time',
                     title='Perturbation', labels=new)
    
    # interpolate the time series, return the data
    mets  = [x for x in sim_data.colnames if x != 'time']
    if time_points:
        interpolated_data = [interp1d(sim_data['time'], sim_data[y], kind='cubic')(time_points) for y in mets]
        df = pd.DataFrame(interpolated_data, index=[model.mapping[x] for x in mets], columns=time_points)
    else:
        time_points = sim_data['time']
        df = pd.DataFrame(sim_data, index=[model.mapping[x] for x in mets], columns=time_points)
    df.columns.name = 'time'
    df.index.name = 'metabolites'
    return df
              
            
def main(n_parameterizations, condition, new_condition='Glucose'):
    """
        
    """
    extended_core = cobra.io.read_sbml_model(settings.ECOLI_EXCORE_FNAME)
    columns = [rxn.id for rxn in extended_core.reactions]
    index = [met.id for met in extended_core.metabolites]
    S = pd.DataFrame(create_stoichiometric_matrix(extended_core), index=index, columns=columns)
    
    ex_rxns = [rxn.id for rxn in extended_core.reactions if rxn.boundary]
    S_sub = remove_mets_rxns(S, ex_rxns)
    
    ## add biomass_e as metabolite, excrete it
    S_sub.loc['biomass_e'] = np.zeros(len(S_sub.columns))
    S_sub.loc['biomass_e', 'Biomass_Ecoli_core_w_GAM'] = 1
    
    ### flux
    flux_solutions = pd.read_pickle(settings.CACHE_DIR+'/flux_solutions.p')
    flux_solution = flux_solutions[condition]
    
    flip=[]
    for rxn in extended_core.reactions:
        if not rxn.boundary:
            if float(flux_solution.loc[rxn.id]) < 0:
                flip += [rxn.id]
                
    for rxn in flip: 
        S_sub.loc[:,rxn] = S_sub.loc[:,rxn] * -1
    flux_solution = abs(flux_solution)
    
    M = ['g6p_c','6pgc_c','r5p_c','s7p_c','fdp_c','dhap_c','g3p_c','pep_c','pyr_c','ac_c',
         'cit_c','akg_c','succ_c','fum_c','mal__L_c','oaa_c','glx_c']
    R = ['PYK','PPS','PYK','PPC','PPCK','ME1','PTAr',
         'PDH','CS','ICDHyr','AKGDH','SUCDi','FUM','MDH','MALS','ICL']
    
    base       = [((),())] 
    single_int = allostery(M,R,1)
    double_int = allostery(M,R,2)
    topologies = base + single_int + double_int

    ### time points to save and number of parameterizations per model
    time_points = [0, 10, 15, 25, 40, 45, 55, 70]

    job_index = int(os.environ.get("LSB_JOBINDEX"))
    number_jobs = int(os.environ.get("LSB_JOBINDEX_END"))
    chunk = topologies[job_index-1::number_jobs]
        
    ### iterate over the topologies in this chunk
    results = dict()
    for j in range(0,len(chunk)):
        model = MakeModel(S_sub.copy())
        for rxn in model.rates:
            model.add_rate(rxn)
        model.add_regulation(chunk[j])
        
        ### iterate over the parameterizations
        result = list()
        for _ in range(n_parameterizations): # multiple parameterizations per model topology
            model.parameterize(flux_solution)
            model.to_antimony()
            
            model = perturb(model, condition, new_condition)
            df_sim = simulate_antimony_model(model, plot=True, time_points = time_points)
            result += [df_sim, pd.Series(model.P)]
            
        results[model.tag] = result

    with open(settings.SERVER_SCRATCH+'/part'+str(job_index)+'.p', 'wb') as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL) 

if __name__ == '__main__':
    #condition = 'Succinate' # initial carbon source, switching to glucose
    #n_parameterizations = 1
    condition = str(sys.argv[1]) 
    n_parameterizations = int(sys.argv[2])
    main(n_parameterizations, condition)

# test: bsub -J "array[1-200]" -W 120:00 -R "rusage[mem=5000]" python ensemble.py Pyruvate 2
# delta_gibbs = Delta_gibbs(settings.DATA_DIR+'/Gerosa_2015_S2.xlsx')

