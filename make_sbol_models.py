from typing import Tuple

import sbol3
import tyto

from helpers import *

MODEL_FILE = '210919_try01.nt'
PROJECT_NAMESPACE = 'http://bbn.com/crispr-kill-switch'


###########################
# Modules:

def constitutive(target: sbol3.Feature) -> sbol3.Feature:
    """Add a constitutive promoter regulating the target feature

    :param target: CDS or ncRNA to regulate
    :return: newly created constitutive promoter
    """
    system: sbol3.Component = get_toplevel(target)
    containers = [c.subject for c in system.constraints
                  if c.restriction == sbol3.SBOL_CONTAINS and c.object == target.identity]
    if len(containers) != 1:
        raise ValueError(f'Should be precisely one container of a constitutive target, but found {len(containers)}')
    vector = containers[0].lookup()
    promoter = add_subfeature(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.promoter]))
    regulate(promoter, target)
    return promoter


def make_crispr_module(vector: sbol3.Feature) -> sbol3.Feature:
    """Add a CRISPR module to the system, comprising both genome editing and kill switch

    :param vector: Vector into which the coding materials for the CRISPR module will be added
    :return: tuple of sgRNA1 and sgRNA2 ncRNA regions, for attaching regulation to
    """
    # find system containing the vector
    system: sbol3.Component = get_toplevel(vector)

    # Add constitutive Cas9 expression # TODO: Change so that it isn't always constitutive
    cas9_cds = add_subfeature(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.CDS], name="Cas9-coding"))
    constitutive(cas9_cds)
    cas9 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name="Cas9"))
    add_interaction(system, sbol3.SBO_GENETIC_PRODUCTION, {cas9_cds: sbol3.SBO_TEMPLATE, cas9: sbol3.SBO_PRODUCT})

    # Add the sgRNA coding regions
    sgRNA1_dna = add_subfeature(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.sgRNA], name="sgRNA1-coding"))
    sgRNA2_dna = add_subfeature(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.sgRNA], name="sgRNA2-coding"))
    constitutive(sgRNA2_dna)

    # Then their products and binding to Cas9
    sgRNA1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_RNA], name="sgRNA1"))
    sgRNA2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_RNA], name="sgRNA2"))
    add_interaction(system, sbol3.SBO_GENETIC_PRODUCTION, {sgRNA1_dna: sbol3.SBO_TEMPLATE, sgRNA1: sbol3.SBO_PRODUCT})
    add_interaction(system, sbol3.SBO_GENETIC_PRODUCTION, {sgRNA2_dna: sbol3.SBO_TEMPLATE, sgRNA2: sbol3.SBO_PRODUCT})
    Cas9_sgRNA1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_NON_COVALENT_COMPLEX], name="Cas9-sgRNA1"))
    Cas9_sgRNA2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_NON_COVALENT_COMPLEX], name="Cas9-sgRNA2"))
    add_interaction(system, sbol3.SBO_NON_COVALENT_BINDING, name='Cas-gRNA binding',
                    participants={sgRNA1: sbol3.SBO_REACTANT, cas9: sbol3.SBO_REACTANT, Cas9_sgRNA1: sbol3.SBO_PRODUCT})
    add_interaction(system, sbol3.SBO_NON_COVALENT_BINDING, name='Cas-gRNA binding',
                    participants={sgRNA2: sbol3.SBO_REACTANT, cas9: sbol3.SBO_REACTANT, Cas9_sgRNA2: sbol3.SBO_PRODUCT})

    # Finally, the Cas9 complex editing actions
    genome = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='genome'))
    add_interaction(system, tyto.SBO.cleavage, name='Cas cleavage',
                    participants={Cas9_sgRNA1: sbol3.SBO_REACTANT, vector: sbol3.SBO_REACTANT})
    add_interaction(system, tyto.SBO.cleavage, name='Cas cleavage',
                    participants={Cas9_sgRNA2: sbol3.SBO_REACTANT, genome: sbol3.SBO_REACTANT})
    add_interaction(system, sbol3.SBO_DEGRADATION, name='Cas degradation',
                    participants={Cas9_sgRNA1: sbol3.SBO_REACTANT})
    add_interaction(system, sbol3.SBO_DEGRADATION, name='Cas degradation',
                    participants={Cas9_sgRNA2: sbol3.SBO_REACTANT})

    # Return the gRNA coding regions for use in establishing regulation
    return sgRNA1_dna


def make_tf_module(system: sbol3.Component, vector: sbol3.Feature, target: sbol3.Feature, repressor: bool): # TODO: Ask Jake why this is a system and not a vector?
    """Add a transcription factor module to the system

    :param system: ???
    :param vector: the AAV genome that this all goes in?
    :param target: thing that is controlled by the TF
    :param repressor: True or false, is this TF a repressor or not
    """

    # I am handed the system

    # Figure out what the vector is
    # TODO: What is the vector that the cds goes into?
    # Hand myself the AAV vector

    # Add constitutive expression of the TF # TODO: Make it not just constitutive
    tf_cds = add_subfeature(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.CDS], name="TF-coding"))
    constitutive(tf_cds)
    tf = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name="TF"))
    add_interaction(system, sbol3.SBO_GENETIC_PRODUCTION, {tf_cds: sbol3.SBO_TEMPLATE, tf: sbol3.SBO_PRODUCT})

    # Make the other promoter
    promoter = add_subfeature(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.promoter]))
    regulate(promoter, target)

    # Add interactions
    # TF activation # TODO: Add if statement to make a repressor option
    add_interaction(system, sbol3.SBO_CONTROL, name='TF Activation',
                    participants={tf: sbol3.SBO_REACTANT, promoter: sbol3.SBO_REACTANT})
    # TF degradation
    add_interaction(system, sbol3.SBO_DEGRADATION, name='TF degradation',
                    participants={tf: sbol3.SBO_REACTANT})

    # Not sure what I am returning
    return promoter

def make_recombinase_module(system: sbol3.Component):
    pass # TODO: implement

###########################
# Actually make the model
doc = sbol3.Document()
sbol3.set_namespace(PROJECT_NAMESPACE)

system = sbol3.Component('Basic_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Basic Kill Switch")
doc.add(system)

aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
sgRNA1_dna = make_crispr_module(aav)
constitutive(sgRNA1_dna)

# Try the TF
system = sbol3.Component('TF_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="TF Kill Switch")
doc.add(system)
sgRNA1_dna = make_crispr_module(aav)
make_tf_module(system, aav, sgRNA1_dna, False)

# Try the Cre
# system = sbol3.Component('Cre_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Cre recombinase Kill Switch")
# doc.add(system)
# sgRNA1_dna = make_recombinase_module(aav)

# Write the model file
doc.write(MODEL_FILE, sbol3.SORTED_NTRIPLES)
