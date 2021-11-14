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

# What happens if we allow 2 TFs or 2 recombinases?
# Chain: A->A, A->R, R->A, R->R, similar for On & Off, for a total of 8
# Joint:
#   Interesting: A+R, On/Off
#   Pointless: A+A, R+R, On/On, Off,Off
# So that would add either 10 or 14 cases (depending on whether we include the pointless ones)
# We'll include the pointless ones (why not?) and consider all 31 cases

# set up combinatorics
TFs = ['Activator', 'Repressor']
recombinases = ['Cre-on', 'Cre-off']
regulator_options = ['Activator', 'Repressor', 'Cre-on', 'Cre-off', None]
orders = ['Chain', 'Joint', '']
combinations = [[r1, r2, order] for r1, r2, order in product(regulator_options, regulator_options, orders)
                if ((len(order) and (r1 and r2)) or (not len(order) and not r2))
                and ((order != 'Joint') or r1 <= r2)]


assert len(combinations) == 31

doc = sbol3.Document()
sbol3.set_namespace(PROJECT_NAMESPACE)

print('Generating models', end='')
for r1, r2, order in combinations:
    # Create system
    name = f'{" ".join(filter(None,[order,r1,r2]))} Kill Switch'.strip()
    system = sbol3.Component(sbol3.string_to_display_id(name), sbol3.SBO_FUNCTIONAL_ENTITY, name=name)
    doc.add(system)

    # Create CRISPR kill switch & regulators
    aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
    sgRNA1_dna, genome = builders.make_crispr_module(aav)
    if r1:
        if r1 in TFs:
            r1_cds, r1_target = builders.make_tf_module(aav, repressor=(r1 == 'Repressor'))
        else:
            r1_cds, r1_target = builders.make_recombinase_module(aav, cre_on=(r1 == 'Cre-on'))
    if r2:
        if r2 in TFs:
            r2_cds, r2_target = builders.make_tf_module(aav, repressor=(r2 == 'Repressor'), second=(r1 in TFs))
        else:
            r2_cds, r2_target = builders.make_recombinase_module(aav, cre_on=(r2 == 'Cre-on'), second=(r1 in recombinases))

    if r2:
        if order == 'Chain':
            constitutive(r1_cds)
            regulate(r1_target, r2_cds)
            regulate(r2_target, sgRNA1_dna)
        elif order == 'Joint':
            constitutive(r1_cds)
            constitutive(r2_cds)
            regulate(r1_target, sgRNA1_dna)
            regulate(r2_target, sgRNA1_dna)
        else:
            raise ValueError(f'Bad order: {order}')
    elif r1:
        constitutive(r1_cds)
        regulate(r1_target, sgRNA1_dna)
    else:
        constitutive(sgRNA1_dna)

    # TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/324
    inputs = [aav, genome]
    if r1 in recombinases:
        inputs.append(r1_target)
    if r2 in recombinases:
        inputs.append(r2_target)
    system.interface = sbol3.Interface(inputs=inputs, outputs=[aav])

    print('.', end='')

print('done')


# Write the model file
doc.write(MODEL_FILE, sbol3.SORTED_NTRIPLES)
