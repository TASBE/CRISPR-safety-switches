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
# Create CRISPR kill switch
aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
sgRNA1_dna, genome = builders.make_crispr_module(aav)
# run it constitutively
builders.constitutive(sgRNA1_dna)
# TODO: Warning will go away after resolution of hhttps://github.com/SynBioDex/pySBOL3/issues/324
system.interface = sbol3.Interface(inputs=[aav, genome], outputs=[aav])

# Try the TF
system = sbol3.Component('TF_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Activator on Kill Switch")
doc.add(system)
# Create CRISPR kill switch
aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
sgRNA1_dna, genome = builders.make_crispr_module(aav)
# regulate with an activator
tf_cds, tf_promoter = builders.make_tf_module(aav, repressor=False)
regulate(tf_promoter, sgRNA1_dna)
builders.constitutive(tf_cds)
# TODO: Warning will go away after resolution of hhttps://github.com/SynBioDex/pySBOL3/issues/324
system.interface = sbol3.Interface(inputs=[aav, genome], outputs=[aav])


# Try the Cre
system = sbol3.Component('Cre_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Cre-activated Kill Switch")
doc.add(system)
# Create CRISPR kill switch
aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
sgRNA1_dna, genome = builders.make_crispr_module(aav)
# regulate with Cre_on
cre_cds, cre_region = builders.make_recombinase_module(aav, cre_on=True)
regulate(cre_region, sgRNA1_dna)
builders.constitutive(cre_cds)
# TODO: Warning will go away after resolution of hhttps://github.com/SynBioDex/pySBOL3/issues/324
system.interface = sbol3.Interface(inputs=[aav, genome, cre_region], outputs=[aav])

# Try chained regulation and dual regulation
system = sbol3.Component('Joint_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Dual-activated Kill Switch")
doc.add(system)
# Create CRISPR kill switch
aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
sgRNA1_dna, genome = builders.make_crispr_module(aav)
# regulate with an activator
tf_cds, tf_promoter = builders.make_tf_module(aav, repressor=False)
regulate(tf_promoter, sgRNA1_dna)
builders.constitutive(tf_cds)
# regulate with Cre_on
cre_cds, cre_region = builders.make_recombinase_module(aav, cre_on=True)
regulate(cre_region, sgRNA1_dna)
builders.constitutive(cre_cds)
# TODO: Warning will go away after resolution of hhttps://github.com/SynBioDex/pySBOL3/issues/324
system.interface = sbol3.Interface(inputs=[aav, genome, cre_region], outputs=[aav])

system = sbol3.Component('Chain_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Repression of Cre-activated Kill Switch")
doc.add(system)
# Create CRISPR kill switch
aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
sgRNA1_dna, genome = builders.make_crispr_module(aav)
# regulate with Cre_on
cre_cds, cre_region = builders.make_recombinase_module(aav, cre_on=True)
regulate(cre_region, sgRNA1_dna)
# regulate Cre-on with a repressor
tf_cds, tf_promoter = builders.make_tf_module(aav, repressor=True)
regulate(tf_promoter, cre_cds)
builders.constitutive(tf_cds)
# TODO: Warning will go away after resolution of hhttps://github.com/SynBioDex/pySBOL3/issues/324
system.interface = sbol3.Interface(inputs=[aav, genome, cre_region], outputs=[aav])

# Eventual full generation should have 17 combinations:
# - 1 unregulated kill switch
# - 4 single regulators kill switches (activator, repressor, Cre_on, Cre_off)
# - 4 combinations of 2 regulators, applied in 3 different ways: TF -> Cre -> gRNA, Cre -> TF -> gRNA, Cre|TF -> gRNA
# Of these we expect half to be viable options, with the other half speeding instead of slowing

# Write the model file
doc.write(SAMPLER_FILE, sbol3.SORTED_NTRIPLES)
