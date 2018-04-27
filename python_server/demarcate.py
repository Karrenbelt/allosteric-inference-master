import settings
import cobra
import re

def main():
    """
    This script reads in the genome scale metabolic model and the core model of 
    E. coli. Then an extended version of the core model is generated, and the 
    namespace of the genome-scale metabolic model is adopted for as far as 
    possible. Additional pathways added or modified: 
        uptake: gluconate, glycerol, galactose, fructose.
        internal: ACCOA synthase, ED pathway, malate oxidase
    """
    GSMM = cobra.io.read_sbml_model(settings.ECOLI_SBML_FNAME)
    core = cobra.io.read_sbml_model(settings.ECOLI_CORE_FNAME)
    outfile = settings.CACHE_DIR+'/extended_core.xml'
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
                
    # change inconsistency in metabolite naming, adopt iJO1366 standard
    manual_rxn_map = {'EX_glc_e':'EX_glc__D_e','EX_gln_L_e':'EX_gln__L_e','EX_glu_L_e':'EX_glu__L_e','EX_lac_D_e':'EX_lac__D_e','EX_mal_L_e':'EX_mal__L_e'}
    for met in core.metabolites:
        if met not in GSMM.metabolites:
            met.id = met.id.replace('_','__',1)
    for rxn in core.reactions:
        if rxn not in GSMM.reactions:
            if rxn.id in manual_rxn_map:
                rxn.id = manual_rxn_map[rxn.id]
    core.repair()
    
    
    # 2. add additional reactions of interest from iJO1366
    
    # add transport reactions in similar style as the others in the core model
    def add_from_iJO1366(rxn_ids_list):
        [core.add_reactions([GSMM.reactions.get_by_id(rxn)]) for rxn in rxn_ids_list]
        return core
    
    ## adding pathways for gluconate
    GLCNt2r = cobra.Reaction('GLCNt2r') # GLCNt2rpp, GLCNtex
    h_e = GSMM.metabolites.get_by_id('h_e')
    h_c = GSMM.metabolites.get_by_id('h_c')
    glcn_e = GSMM.metabolites.get_by_id('glcn_e')
    glcn_c = GSMM.metabolites.get_by_id('glcn_c')
    GLCNt2r.add_metabolites({glcn_e: -1.0, h_e:-1.0, glcn_c: 1.0, h_c:1.0})
    core.add_reactions([GLCNt2r])
    core = add_from_iJO1366(['GNK','EX_glcn_e'])
    
    ## adding pathways for glycerol
    GLYCt = cobra.Reaction('GLYCt') # GLYCtex, GLYCtpp
    glyc_e = GSMM.metabolites.get_by_id('glyc_e')
    glyc_c = GSMM.metabolites.get_by_id('glyc_c')
    GLYCt.add_metabolites({glyc_e: -1.0, glyc_c: 1.0})
    core.add_reactions([GLYCt])
    core = add_from_iJO1366(['GLYCDx','DHAPT','G3PD2','G3PT','GLYK','EX_glyc_e'])
    
    ## adding pathways for galactose
    GALabc = cobra.Reaction('GALabc') # GALtex, GALabspp
    gal_e = GSMM.metabolites.get_by_id('gal_e')
    gal_c = GSMM.metabolites.get_by_id('gal_c')
    atp_c = GSMM.metabolites.get_by_id('atp_c')
    adp_c = GSMM.metabolites.get_by_id('adp_c')
    h2o_c = GSMM.metabolites.get_by_id('h2o_c')
    h_c = GSMM.metabolites.get_by_id('h_c')
    pi_c = GSMM.metabolites.get_by_id('pi_c')
    GALabc.add_metabolites({gal_e: -1.0, atp_c: -1.0, h2o_c: -1.0, gal_c: 1.0, adp_c: 1.0, pi_c: 1.0, h_c: 1.0})
    core.add_reactions([GALabc])
    
    UGLT = cobra.Reaction('UGLT') # effectively removes udpg_c and updgal_c from the original reaction
    gal1p_c = GSMM.metabolites.get_by_id('gal1p_c') 
    g1p_c = GSMM.metabolites.get_by_id('g1p_c')
    UGLT.add_metabolites({gal1p_c: -1.0, g1p_c: 1.0})
    core.add_reactions([UGLT])
    core = add_from_iJO1366(['PGMT','GALKr','EX_gal_e'])
    
    ## changing pathways for fructose
    core.remove_reactions(['FRUpts2']) # change entry point to fpb_c instead of f6p_c (ambiguous), to match Gerosa 2015
    FRUpts = cobra.Reaction('FRUpts')
    fru_e = GSMM.metabolites.get_by_id('fru_e')
    pep_c = GSMM.metabolites.get_by_id('pep_c')
    f1p_c = GSMM.metabolites.get_by_id('f1p_c')
    pyr_c = GSMM.metabolites.get_by_id('pyr_c')
    FRUpts.add_metabolites({fru_e: -1.0, pep_c: -1.0, f1p_c: 1.0, pyr_c: 1.0})
    core.add_reactions([FRUpts])
    core = add_from_iJO1366(['FRUK'])
    
    ## add ED pathway, ACS for acetate, and MOX for malate oxidase (MQO, MDH2, MDH3: all the same reaction, diff cofactor)
    core = add_from_iJO1366(['EDD','EDA'] + ['ACS'] + ['MOX'])
    
    core.reactions.get_by_id('EX_glc__D_e').bounds = [0, 1000] # is not set to zero in the core model initially
    
    core.repair()
    
    ### biomass equation will be changed later
    for rxn in core.reactions:
        if not rxn.boundary:
            if rxn.check_mass_balance():
                print('\nThe reaction following reaction is not mass balanced: {}. Imbalance: {}'.format(rxn.id, rxn.check_mass_balance()) )

    cobra.io.write_sbml_model(core, outfile)                
    print('\n The constraint-based cobra mode has been written, it can be found here: {}\n'.format(outfile) )

if __name__ == '__main__':
    main()
                