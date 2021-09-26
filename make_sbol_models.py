import sbol3
import builders
from helpers import add_feature
from shared_global_names import *

###########################
# Actually make the model
doc = sbol3.Document()
sbol3.set_namespace(PROJECT_NAMESPACE)

system = sbol3.Component('Basic_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Basic Kill Switch")
doc.add(system)

aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
sgRNA1_dna = builders.make_crispr_module(aav)
builders.constitutive(sgRNA1_dna)

# Try the TF
# system = sbol3.Component('TF_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="TF Kill Switch")
# doc.add(system)
# sgRNA1_dna = make_crispr_module(aav)
# make_tf_module(system, aav, sgRNA1_dna, False)

# Try the Cre
# system = sbol3.Component('Cre_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Cre recombinase Kill Switch")
# doc.add(system)
# sgRNA1_dna = make_recombinase_module(aav)

# Write the model file
doc.write(MODEL_FILE, sbol3.SORTED_NTRIPLES)
