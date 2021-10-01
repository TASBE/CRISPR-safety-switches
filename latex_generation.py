import logging
import itertools
from typing import Dict, List, Optional

import sbol3
import tyto
from sbol_utilities.workarounds import id_sort

from sbol_utilities.component import all_in_role, in_role

name_to_symbol = {
    'sgRNA1': '\\gRna{1}',
    'sgRNA2': '\\gRna{2}',
    'Cas9': '\\proSp{\\cas{}}',
    'Cas9-sgRNA1': '\\cplx{\\cas}{1}',
    'Cas9-sgRNA2': '\\cplx{\\cas}{2}',
    'postedit Cas9-sgRNA1': '\\bound{\\cplx{\\cas}{1}}',
    'postedit Cas9-sgRNA2': '\\bound{\\cplx{\\cas}{2}}',
    'genome': '\\hostGen{}',
    'edited genome': '\\edited{\\hostGen{}}',
    'AAV': '\\vectorGen{}',
    'Cas-gRNA binding': '\\gRnaBind{}',
    'Cas degradation': '\\casCompDegradeRate{}',
    'Cas cleavage': '\\casCutRate{}',
    'TF': '\\proSp{TF}'
}
"""Dictionary mapping from SBOL names to LaTeX symbols in our convention"""


def maybe_concentration(feature: sbol3.Feature) -> str:
    """Determine whether we are working with a concentration or a count based on type

    :param feature: Feature to be evaluated
    :return: symbol, possibly wrapped in a concentration
    """
    if feature.name not in name_to_symbol:
        raise ValueError(f'No symbol known for name: "{feature.name}"')
    symbol = name_to_symbol[feature.name]
    if sbol3.SBO_DNA in feature.types:
        return symbol
    else:
        return f'\\conc{{{symbol}}}'


def differential(feature: sbol3.Feature) -> str:
    """Return the "dX/dt" term of the ODE and equals sign

    :param feature: feature to get a differential for
    :return: LaTeX string
    """
    return f'\\diff{{{maybe_concentration(feature)}}}{{t}} & = '


def regulation_term(interaction: sbol3.Interaction) -> str:
    """Generate a term for regulation by transcription factor or recombinase

    :param interaction: Regulation interaction to serialize
    :return: LaTeX serialization
    """
    # Need i_type to see what type of regulation is happening
    i_type = interaction.types[0]
    # Make TF Equations
    if i_type == sbol3.SBO_INHIBITION:
        regulator_species = name_to_symbol['TF']
        # TODO: Replace K and n with variables
        return f'\\frac{{(K_R)^n}}{{(K_R)^n + {regulator_species}^n}}' # TODO: Use maybe_concentrations to add the [], instead of hardcode?
        return(regulation_term)
    elif i_type == sbol3.SBO_STIMULATION:
        regulator_species = name_to_symbol['TF']
        # TODO: Replace K and n with variables
        return f'\\frac{{[{regulator_species}]^n}}{{(K_A)^n + [{regulator_species}]^n}}'
    # Make Cre equations
    # elif i_type == cre_recombinase: # TODO: Implement Cre
    #     pass
    else:
        logging.warning(f'Cannot serialize interaction {interaction.identity} of type {tyto.SBO.get_term_by_uri(i_type)}')
        return ''


def interaction_to_term(feature: sbol3.Feature, interaction: sbol3.Interaction,
                        regulation: Dict[sbol3.Feature, List[sbol3.Interaction]],
                        containers: Dict[sbol3.Feature, List[sbol3.Feature]]) -> Optional[str]:
    """Generate an equation term for a given interaction, with respect to the included feature

    :param feature: Target of the term
    :param interaction: Interaction to get an equation for
    :param regulation: Dictionary of regulation interactions in the system
    :param containers: Dictionary of container relationships in system
    :return: LaTeX equation term
    """
    if len(interaction.types) != 1:
        raise ValueError(f'Expected 1 interaction type but found {len(interaction.types)} in {interaction.identity}')
    if len(feature.types) != 1:
        raise ValueError(f'Expected 1 feature type but found {len(feature.types)} in {feature.identity}')
    # find the participation for this feature and its role therein
    feature_participation = [p for p in interaction.participations if p.participant == feature.identity]
    if len(feature_participation) != 1:
        raise ValueError(f'Expected feature in 1 participant, but found {len(feature_participation)} in {interaction.identity}')
    if len(feature_participation[0].roles) != 1:
        raise ValueError(f'Do not know how to serialize multi-role participation {feature_participation[0]}')
    i_type = interaction.types[0]
    f_type = feature.types[0]
    role = feature_participation[0].roles[0]
    # serialized based on interaction type and role
    if i_type == sbol3.SBO_GENETIC_PRODUCTION:
        if role == sbol3.SBO_TEMPLATE:
            return None  # templates don't get equations - they are taken as regulator for products
        elif role == sbol3.SBO_PRODUCT:
            species = name_to_symbol[feature.name]
            template = in_role(interaction,sbol3.SBO_TEMPLATE)
            # modulation is the regulation of either the template or the product
            modulation = ''.join(regulation_term(r) for r in regulation[feature] + regulation[template])
            # context is the constraints of the template
            context = ''.join(maybe_concentration(ct) for ct in containers[template])
            if f_type == sbol3.SBO_RNA:
                prod_rate = f'\\txRate{{{species}}}'
                deg_rate = f'\\rnaDegradeRate{{}}'
            elif f_type == sbol3.SBO_PROTEIN:
                prod_rate = f'\\txtlRate{{{species}}}'
                deg_rate = f'\\proDegradeRate{{{species}}}'
            else:
                raise ValueError(f'Cannot handle type {tyto.SBO.get_term_by_uri(f_type)} in {feature_participation[0]}')
            return f'+ {prod_rate}{modulation}{context} - {deg_rate}{maybe_concentration(feature)}'
        else:
            logging.warning(f'Cannot serialize role in {interaction.identity} of type {tyto.SBO.get_term_by_uri(i_type)}')
    elif i_type == tyto.SBO.cleavage:
        if interaction.name == 'Cas cleavage':
            reactants = [maybe_concentration(f) for f in all_in_role(interaction, sbol3.SBO_REACTANT)]
            rate = '\\casCutRate{{}}'
            if role == sbol3.SBO_REACTANT:
                sign = '-'
            elif role == sbol3.SBO_PRODUCT:
                sign = '+'
            else:
                raise ValueError (f'Unexpected role in {interaction.identity}: {tyto.SBO.get_term_by_uri(role)}')
            return f'{sign} {rate}' + ''.join(reactants)
        else:
            raise ValueError(f'No model for cleavage {interaction.name} in {interaction.identity}')
    elif i_type == sbol3.SBO_DEGRADATION:
        if len(interaction.participations) != 1:
            raise ValueError(f'Degradation assumed to have 1 participant, found {len(interaction.participations)} in {interaction.identity}')
        if f_type == sbol3.SBO_RNA:
            deg_rate = f'\\rnaDegradeRate{{}}'
        elif f_type == sbol3.SBO_PROTEIN:
            species = name_to_symbol[feature.name]
            deg_rate = f'\\proDegradeRate{{{species}}}'
        else:
            deg_rate = name_to_symbol[interaction.name]

        return f'- {deg_rate}{maybe_concentration(feature)}'

        pass
    elif i_type == sbol3.SBO_NON_COVALENT_BINDING:
        reactants = [maybe_concentration(f) for f in all_in_role(interaction, sbol3.SBO_REACTANT)]
        rate = name_to_symbol[interaction.name] # TODO: move this into actual parameters rather than name
        if role == sbol3.SBO_REACTANT:
            sign = '-'
        elif role == sbol3.SBO_PRODUCT:
            sign = '+'
        else:
            raise ValueError(f'Cannot handle type {tyto.SBO.get_term_by_uri(f_type)} in {interaction.identity}')
        return f'{sign} {rate}' + ''.join(reactants)
    elif i_type == sbol3.SBO_INHIBITION or i_type == sbol3.SBO_STIMULATION:
        # Pass for the regulation interactions that are taken care of in the other function, so you don't get a warning
        pass
    else:
        logging.warning(f'Cannot serialize interaction {interaction.identity} of type {tyto.SBO.get_term_by_uri(i_type)}')
        return None


def make_latex_model(system: sbol3.Component) -> str:
    """Generate a set of LaTeX equations for the identified system:

    :param system: system for which a model is to be generated
    :return: string serialization of LaTeX equation collection
    """
    # for each feature, collect all of the interactions and constraints that it participates in
    interactions = {f: [i for i in system.interactions
                           if [p for p in i.participations if p.participant == f.identity]]
                       for f in system.features}
    containers = {f: [c.subject.lookup() for c in system.constraints if c.restriction == sbol3.SBOL_CONTAINS and c.object == f.identity]
                  for f in system.features}
    regulators = {f: [c.subject.lookup() for c in system.constraints if c.restriction == sbol3.SBOL_MEETS and c.object == f.identity]
                  for f in system.features}
    regulation = {f: list(itertools.chain(*(interactions[r] for r in regulators[f]))) for f in regulators}

    # generate an ODE based on the roles in the interactions
    equation_latex = []
    for f in id_sort(system.features):
        interaction_terms = [t for t in [interaction_to_term(f, i, regulation, containers) for i in id_sort(interactions[f])] if t]
        # If there is at least one term, then add an equation
        if interaction_terms:
            diff_term = differential(f)
            equation_latex.append(diff_term + " ".join(interaction_terms).removeprefix("+"))

    ## Generate the actual document
    # write section header
    latex =  f'\\subsection{{{system.name or system.display_id}}}\n\\label{{s:{system.display_id}}}\n'
    latex += f'% Equations generated from {system.identity}\n\n'
    if system.description:
        latex += f'{system.description}\n\n'
    # write equations
    latex += f'\\begin{{align}}\n'
    latex += '\\\\\n'.join(equation_latex)
    latex += f'\n\\end{{align}}\n\n'

    return latex
