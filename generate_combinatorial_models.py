from itertools import product

import sbol3

import builders
from sbol_utilities.component import add_feature, regulate, constitutive
from shared_global_names import *

# Full generation should have 17 combinations:
# - 1 unregulated kill switch
# - 4 single regulators kill switches (activator, repressor, Cre_on, Cre_off)
# - 4 combinations of 2 regulators, applied in 3 different ways: TF -> Cre -> gRNA, Cre -> TF -> gRNA, Cre|TF -> gRNA
# Of these we expect half to be viable options, with the other half speeding instead of slowing

# set up combinatorics
TF_options = ['Activator', 'Repressor', None]
Cre_options = ['Cre-on', 'Cre-off', None]
orders = ['Chain-Cre-TF', 'Chain-TF-Cre', 'Joint', '']
combinations = [[tf, cre, order] for tf, cre, order in product(TF_options, Cre_options, orders)
                if (len(order) and (tf and cre)) or (not len(order) and not (tf and cre))]


assert len(combinations) == 17

doc = sbol3.Document()
sbol3.set_namespace(PROJECT_NAMESPACE)

print('Generating models', end='')
for tf, cre, order in combinations:
    # Create system
    name = f'{" ".join(filter(None,[order,tf,cre]))} Kill Switch'.strip()
    system = sbol3.Component(sbol3.string_to_display_id(name), sbol3.SBO_FUNCTIONAL_ENTITY, name=name)
    doc.add(system)

    # Create CRISPR kill switch & regulators
    aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
    sgRNA1_dna, genome = builders.make_crispr_module(aav)
    if tf:
        tf_cds, tf_promoter = builders.make_tf_module(aav, repressor=(tf == 'Repressor'))
    if cre:
        cre_cds, cre_region = builders.make_recombinase_module(aav, cre_on=(cre == 'Cre-on'))

    if tf and cre:
        if order == 'Chain-Cre-TF':
            constitutive(cre_cds)
            regulate(cre_region, tf_cds)
            regulate(tf_promoter, sgRNA1_dna)
        elif order == 'Chain-TF-Cre':
            constitutive(tf_cds)
            regulate(tf_promoter, cre_cds)
            regulate(cre_region, sgRNA1_dna)
        elif order == 'Joint':
            constitutive(cre_cds)
            constitutive(tf_cds)
            regulate(tf_promoter, sgRNA1_dna)
            regulate(cre_region, sgRNA1_dna)
        else:
            raise ValueError(f'Bad order: {order}')
    elif tf:
        constitutive(tf_cds)
        regulate(tf_promoter, sgRNA1_dna)
    elif cre:
        constitutive(cre_cds)
        regulate(cre_region, sgRNA1_dna)
    else:
        constitutive(sgRNA1_dna)

    # TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/315
    # TODO: interfaces will change to interface after resolution of https://github.com/SynBioDex/pySBOL3/issues/316
    if cre:
        system.interfaces = sbol3.Interface(input=[aav, genome, cre_region], output=[aav])
    else:
        system.interfaces = sbol3.Interface(input=[aav, genome], output=[aav])

    print('.', end='')

print('done')


# Write the model file
doc.write(MODEL_FILE, sbol3.SORTED_NTRIPLES)
