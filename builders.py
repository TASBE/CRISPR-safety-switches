from typing import Tuple

import sbol3
import tyto

###########################
# Modules:
from sbol_utilities.workarounds import get_toplevel
from sbol_utilities.component import constitutive, contains, add_feature, add_interaction, order

from shared_global_names import RECOMBINATION


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
    add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={cas9: sbol3.SBO_REACTANT})

    # Add the sgRNA coding regions
    sgRNA1_dna = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.sgRNA], name="sgRNA1-coding"))
    sgRNA2_dna = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.sgRNA], name="sgRNA2-coding"))
    constitutive(sgRNA2_dna)

    # Then their products and binding to Cas9
    sgRNA1 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_RNA], name="sgRNA1"))
    sgRNA2 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_RNA], name="sgRNA2"))
    add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {sgRNA1_dna: sbol3.SBO_TEMPLATE, sgRNA1: sbol3.SBO_PRODUCT})
    add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {sgRNA2_dna: sbol3.SBO_TEMPLATE, sgRNA2: sbol3.SBO_PRODUCT})
    add_interaction(sbol3.SBO_DEGRADATION, name='gRNA degradation', participants={sgRNA1: sbol3.SBO_REACTANT})
    add_interaction(sbol3.SBO_DEGRADATION, name='gRNA degradation', participants={sgRNA2: sbol3.SBO_REACTANT})
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


def make_tf_module(vector: sbol3.Feature, repressor: bool, second: bool = False) -> Tuple[sbol3.Feature, sbol3.Feature]:
    """Add a transcription factor regulation module to the system

    :param vector: Vector into which the coding materials for the TF module will be added
    :param repressor: true for repressor, false for activator
    :param second: true if this is the second, and thus should be "TF2" instead of "TF"
    :returns: tuple of CDS and promoter features, for connecting to regulation
    """

    # find system containing the vector
    system = get_toplevel(vector)
    if not isinstance(system, sbol3.Component):
        raise ValueError(f'System should be a component but was not: {system}')
    name = "TF2" if second else "TF"

    # Add the cds of the TF, the TF, and the production relation between them
    tf_cds = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.CDS], name=f'{name}-coding'))
    tf = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name=name))
    add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {tf_cds: sbol3.SBO_TEMPLATE, tf: sbol3.SBO_PRODUCT})
    add_interaction(sbol3.SBO_DEGRADATION, name=f'{name} degradation', participants={tf: sbol3.SBO_REACTANT})

    # Make the promoter that is regulated by the TF and add its regulation
    promoter = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.promoter]))
    if repressor:
        add_interaction(sbol3.SBO_INHIBITION, name=f'TF Repression',
                        participants={tf: sbol3.SBO_INHIBITOR, promoter: sbol3.SBO_INHIBITED})
    else:
        add_interaction(sbol3.SBO_STIMULATION, name=f'TF Activation',
                        participants={tf: sbol3.SBO_STIMULATOR, promoter: sbol3.SBO_STIMULATED})

    # Return the cds and the promoter
    return tf_cds, promoter


def make_recombinase_module(vector: sbol3.Feature, cre_on: bool, second: bool = False) -> Tuple[sbol3.Feature, sbol3.Feature]:
    """Add a Cre-recombinase regulation module to the system

    :param vector: Vector into which the coding materials for the Cre module will be added
    :param cre_on: true if Cre activates expression, false if Cre shuts off expression
    :param second: true if this is the second, and thus should be "CreH" instead of "Cre"
    :returns: tuple of CDS and regulatory features, for connecting to regulation
    """

    # find system containing the vector
    system = get_toplevel(vector)
    if not isinstance(system, sbol3.Component):
        raise ValueError(f'System should be a component but was not: {system}')
    name = "CreH" if second else "Cre"

    # Add the cds of the TF, the TF, and the production relation between them
    cre_cds = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.CDS], name=f'{name}-coding'))
    cre = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name=name))
    add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {cre_cds: sbol3.SBO_TEMPLATE, cre: sbol3.SBO_PRODUCT})
    add_interaction(sbol3.SBO_DEGRADATION, name=f'{name} degradation', participants={cre: sbol3.SBO_REACTANT})

    # Make the promoter region that is regulated by the TF and add its regulation
    cre_region = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.engineered_region], name=f'{name} regulated region'))
    edited_cre_region = contains(vector, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.engineered_region], name=f'edited {name} regulated region'))
    promoter = contains(cre_region, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.promoter], name=f'{name} region promoter'))
    cre_target1 = contains(cre_region, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.binding_site], name=f'{name} 5\' target'))
    cre_target2 = contains(cre_region, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[tyto.SO.binding_site], name=f'{name} 3\' target'))
    if cre_on:  # wrap the targets around a terminator for Cre to turn off expression
        terminator = contains(cre_region, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[sbol3.SO_TERMINATOR], name=f'{name}-targeted terminator'))
        order(promoter, cre_target1)
        order(cre_target1, terminator)
        order(terminator, cre_target2)
        add_interaction(RECOMBINATION, name=f'Cre recombination',
                        participants={cre: sbol3.SBO_MODIFIER, cre_region: sbol3.SBO_REACTANT,
                                      terminator: sbol3.SBO_MODIFIED, edited_cre_region: sbol3.SBO_PRODUCT})
    else:  # wrap the targets around the promoter for Cre to turn off expression
        order(cre_target1, promoter)
        order(promoter, cre_target2)
        add_interaction(RECOMBINATION, name=f'Cre recombination',
                        participants={cre: sbol3.SBO_MODIFIER, cre_region: sbol3.SBO_REACTANT,
                                      promoter: sbol3.SBO_MODIFIED, edited_cre_region: sbol3.SBO_PRODUCT})

    # Return the cds and the promoter
    return cre_cds, cre_region
