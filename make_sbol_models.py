import sbol3
import builders
from sbol_utilities.component import add_feature, regulate
from shared_global_names import *

###########################
# Actually make the model
doc = sbol3.Document()
sbol3.set_namespace(PROJECT_NAMESPACE)

# CRISPR only
# system = sbol3.Component('Basic_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Basic Kill Switch")
# doc.add(system)
# aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
# sgRNA1_dna, sgRNA1_rna = builders.make_crispr_module(aav)
# builders.constitutive(sgRNA1_dna)

system = sbol3.Component('Basic_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Basic Kill Switch")
doc.add(system)

aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
sgRNA1_dna, genome = builders.make_crispr_module(aav)
builders.constitutive(sgRNA1_dna)

# TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/315
interface = sbol3.Interface(input=[aav, genome], output=[aav])
# TODO: interfaces will change to interface after resolution of https://github.com/SynBioDex/pySBOL3/issues/316
system.interfaces = interface

# Try the TF
system = sbol3.Component('TF_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Repressor on Kill Switch")
doc.add(system)
aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
sgRNA1_dna, genome = builders.make_crispr_module(aav)
tf_cds, tf_promoter = builders.make_tf_module(aav, False)
regulate(tf_promoter, sgRNA1_dna)
builders.constitutive(tf_cds)

# Try the Cre
# system = sbol3.Component('Cre_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Cre recombinase Kill Switch")
# doc.add(system)
# sgRNA1_dna = make_recombinase_module(aav)

# Write the model file
doc.write(MODEL_FILE, sbol3.SORTED_NTRIPLES)
