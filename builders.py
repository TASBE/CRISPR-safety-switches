from typing import Tuple

import sbol3
import tyto

###########################
# Modules:
from sbol_utilities.workarounds import get_toplevel
from sbol_utilities.component import constitutive, contains, add_feature, add_interaction


def make_crispr_module(vector: sbol3.Feature) -> Tuple[sbol3.Feature, sbol3.Feature]:
    """Add a CRISPR module to the system, comprising both genome editing and kill switch

    :param vector: Vector into which the coding materials for the CRISPR module will be added
    :return: tuple of sgRNA1 ncRNA, genome for attaching regulation to
    """
    # find system containing the vector
    system = get_toplevel(vector)
    if not isinstance(system, sbol3.Component):
        raise ValueError(f'System should be a component but was not: {system}')

    # Add constitutive Cas9 expression # TODO: Change so that it isn't always constitutive
    cas9_cds = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.CDS], name="Cas9-coding"))
    constitutive(cas9_cds)
    cas9 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name="Cas9"))
    add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {cas9_cds: sbol3.SBO_TEMPLATE, cas9: sbol3.SBO_PRODUCT})

    # Add the sgRNA coding regions
    sgRNA1_dna = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.sgRNA], name="sgRNA1-coding"))
    sgRNA2_dna = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.sgRNA], name="sgRNA2-coding"))
    constitutive(sgRNA2_dna)

    # Then their products and binding to Cas9
    sgRNA1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_RNA], name="sgRNA1"))
    sgRNA2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_RNA], name="sgRNA2"))
    add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {sgRNA1_dna: sbol3.SBO_TEMPLATE, sgRNA1: sbol3.SBO_PRODUCT})
    add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {sgRNA2_dna: sbol3.SBO_TEMPLATE, sgRNA2: sbol3.SBO_PRODUCT})
    Cas9_sgRNA1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_NON_COVALENT_COMPLEX], name="Cas9-sgRNA1"))
    Cas9_sgRNA2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_NON_COVALENT_COMPLEX], name="Cas9-sgRNA2"))
    add_interaction(sbol3.SBO_NON_COVALENT_BINDING, name='Cas-gRNA binding',
                    participants={sgRNA1: sbol3.SBO_REACTANT, cas9: sbol3.SBO_REACTANT, Cas9_sgRNA1: sbol3.SBO_PRODUCT})
    add_interaction(sbol3.SBO_NON_COVALENT_BINDING, name='Cas-gRNA binding',
                    participants={sgRNA2: sbol3.SBO_REACTANT, cas9: sbol3.SBO_REACTANT, Cas9_sgRNA2: sbol3.SBO_PRODUCT})

    # Finally, the Cas9 complex editing actions, including "expended" post-edit Cas9
    genome = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='genome'))
    ex_Cas9_1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_NON_COVALENT_COMPLEX], name="postedit Cas9-sgRNA1"))
    ex_Cas9_2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_NON_COVALENT_COMPLEX], name="postedit Cas9-sgRNA2"))
    edited_genome = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='edited genome'))
    add_interaction(tyto.SBO.cleavage, name='Cas cleavage',
                    participants={Cas9_sgRNA1: sbol3.SBO_REACTANT, vector: sbol3.SBO_REACTANT,
                                  ex_Cas9_1: sbol3.SBO_PRODUCT})
    add_interaction(tyto.SBO.cleavage, name='Cas cleavage',
                    participants={Cas9_sgRNA2: sbol3.SBO_REACTANT, genome: sbol3.SBO_REACTANT,
                                  edited_genome: sbol3.SBO_PRODUCT, ex_Cas9_2: sbol3.SBO_PRODUCT})
    add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={Cas9_sgRNA1: sbol3.SBO_REACTANT})
    add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={Cas9_sgRNA2: sbol3.SBO_REACTANT})
    add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={ex_Cas9_1: sbol3.SBO_REACTANT})
    add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={ex_Cas9_2: sbol3.SBO_REACTANT})

    # Return the kill-switch gRNA coding region for use in establishing regulation, genome for output
    return sgRNA1_dna, genome


def make_tf_module(vector: sbol3.Feature, repressor: bool) -> Tuple[sbol3.Feature, sbol3.Feature]:
    """Add a transcription factor regulation module to the system

    :param vector: Vector into which the coding materials for the TF module will be added
    :param repressor: true for repressor, false for activator
    :returns: tuple of CDS and promoter features, for connecting to regulation
    """

    # find system containing the vector
    system = get_toplevel(vector)
    if not isinstance(system, sbol3.Component):
        raise ValueError(f'System should be a component but was not: {system}')

    # Add the cds of the TF, the TF, and the production relation between them
    tf_cds = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.CDS], name="TF-coding"))
    tf = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name="TF"))
    add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {tf_cds: sbol3.SBO_TEMPLATE, tf: sbol3.SBO_PRODUCT})

    # Make the promoter that is regulated by the TF and add its regulation
    promoter = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.promoter]))
    if repressor:
        add_interaction(sbol3.SBO_INHIBITION, name='TF Repression',
                        participants={tf: sbol3.SBO_INHIBITOR, promoter: sbol3.SBO_INHIBITED})
    else:
        add_interaction(sbol3.SBO_STIMULATION, name='TF Activation',
                        participants={tf: sbol3.SBO_STIMULATOR, promoter: sbol3.SBO_STIMULATED})

    # Return the cds and the promoter
    return tf_cds, promoter


def make_recombinase_module(system: sbol3.Component):
    pass  # TODO: implement
