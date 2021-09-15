from typing import Tuple

import sbol3
import tyto

from helpers import *

MODEL_FILE = 'kill_switch_models.nt'
PROJECT_NAMESPACE = 'http://bbn.com/crispr-kill-switch'


###########################
# Modules:

def constitutive(target: sbol3.Feature) -> sbol3.Feature:
    """Add a constitutive promoter regulating the target feature

    :param target: CDS or ncRNA to regular
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

    # Add constitutive Cas9 expression
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


def make_tf_module(system: sbol3.Component, repressor: bool):
    pass # TODO: implement


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

doc.write(MODEL_FILE, sbol3.SORTED_NTRIPLES)
