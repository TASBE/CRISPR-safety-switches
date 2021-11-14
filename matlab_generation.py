import logging
import itertools
from collections import UserDict
from typing import Dict, List, Optional, Union, Tuple

import sbol3
import tyto
from sbol_utilities.helper_functions import id_sort
from sbol_utilities.component import in_role, all_in_role

from shared_global_names import RECOMBINATION


class VariableDictionary(UserDict):
    """Collection of variables, as a wrapper around a dictionary that auto-adds names for missing keys"""

    def __getitem__(self, key):
        if key not in self.data:
            self.data[key] = matlab_name(key)
        return self.data[key]


class ParameterDictionary(UserDict):
    """Collection of variables, as a wrapper around a dictionary that auto-adds names for missing keys"""

    def __getitem__(self, key):
        if key not in self.data:
            if isinstance(key, str):
                self.data[key] = sbol3.string_to_display_id(key)
            else:
                self.data[key] = matlab_name(key)  # TODO: in matlab_name, Interactions are not Features; fix here or there
        return self.data[key]


def matlab_name(feature: sbol3.Feature) -> str:
    """Get a Matlab-compatible variable name for an object

    :param feature: feature to get a Matlab variable name
    :return: Matlab string
    """
    return sbol3.string_to_display_id(feature.name)


def differential(variable: Union[sbol3.Feature, str]) -> str:
    """Return the "dX/dt" variable of the ODE

    :param variable: feature to get a differential for
    :return: Matlab string
    """
    if isinstance(variable, sbol3.Feature):
        variable = matlab_name(variable)
    return f'd_{variable}'


def regulation_term(interaction: sbol3.Interaction, parameters: ParameterDictionary,
                    variables: VariableDictionary) -> str:
    """Generate a term for regulation by transcription factor or recombinase

    :param interaction: Regulation interaction to serialize
    :param parameters: Known parameters for system
    :param variables: Known variables for system
    :return: Matlab equation term
    """
    # Need i_type to see what type of regulation is happening
    i_type = interaction.types[0]
    # Make TF Equations
    if i_type == sbol3.SBO_INHIBITION:
        species = variables[in_role(interaction, sbol3.SBO_INHIBITOR)]
        # TODO: Consider replacing K and n with variables
        k = parameters['K_R']
        n = parameters['n']
        return f'({k}^{n})/({k}^{n} + {species}^{n})'
    elif i_type == sbol3.SBO_STIMULATION:
        species = variables[in_role(interaction, sbol3.SBO_STIMULATOR)]
        # TODO: Consider replacing K and n with variables
        k = parameters['K_A']
        n = parameters['n']
        return f'({species}^{n})/({k}^{n} + {species}^{n})'
    # Make Cre equations
    elif i_type == RECOMBINATION:
        target = in_role(interaction, sbol3.SBO_MODIFIED)
        if any(tyto.SO.promoter.is_ancestor_of(r) for r in target.roles):  # Cre-off
            original = in_role(interaction, sbol3.SBO_REACTANT)
            return f'({variables[original]}/AAV)' # TODO: replace AAV w. variable
        elif any(tyto.SO.terminator.is_ancestor_of(r) for r in target.roles):  # Cre-on
            recombined = in_role(interaction, sbol3.SBO_PRODUCT)
            return f'({variables[recombined]}/AAV)' # TODO: replace AAV w. variable
        else:
            raise ValueError(f'Cannot give term for recombination on roles {target.roles} in {interaction.identity}')
    else:
        logging.warning(f'Cannot serialize regulation {interaction.identity}, type {tyto.SBO.get_term_by_uri(i_type)}')
        return ''


def interaction_to_term(feature: sbol3.Feature, interaction: sbol3.Interaction,
                        regulation: Dict[sbol3.Feature, List[sbol3.Interaction]],
                        containers: Dict[sbol3.Feature, List[sbol3.Feature]], parameters: ParameterDictionary,
                        variables: VariableDictionary) -> Optional[str]:
    """Generate an equation term for a given interaction, with respect to the included feature

    :param feature: Target of the term
    :param interaction: Interaction to get an equation for
    :param regulation: Dictionary of regulation interactions in system
    :param containers: Dictionary of container relationships in system
    :param parameters: Known parameters for system
    :param variables: Known variables for system
    :return: Matlab equation term
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

    # serialize based on interaction type and role
    if i_type == sbol3.SBO_GENETIC_PRODUCTION:
        if role == sbol3.SBO_TEMPLATE:
            return None  # templates don't get equations
        elif role == sbol3.SBO_PRODUCT:
            species = variables[feature]
            template = in_role(interaction, sbol3.SBO_TEMPLATE)
            # modulation is the regulation of either the template or the product
            modulation = '*'.join(regulation_term(r, parameters, variables)
                                 for r in regulation[feature] + regulation[template])
            # context is the constraints of the template
            context = ''.join(variables[ct] for ct in containers[template])
            if f_type == sbol3.SBO_RNA:
                prod_rate = parameters[f'alpha_r_{species}']
            elif f_type == sbol3.SBO_PROTEIN:
                prod_rate = parameters[f'alpha_p_{species}']
            else:
                raise ValueError(f'Cannot handle type {tyto.SBO.get_term_by_uri(f_type)} in {feature_participation[0]}')
            return f'+ {"*".join(filter(None, [prod_rate, modulation, context]))}'
        else:
            logging.warning(f'Cannot serialize role in {interaction.identity}, type {tyto.SBO.get_term_by_uri(i_type)}')
    elif i_type == tyto.SBO.cleavage:
        if interaction.name == 'Cas cleavage':
            reactants = [variables[f] for f in all_in_role(interaction, sbol3.SBO_REACTANT)]
            [variables[f] for f in all_in_role(interaction, sbol3.SBO_PRODUCT)] # Get products into the variable table
            rate = parameters['k_cat']
            if role == sbol3.SBO_REACTANT:
                sign = '-'
            elif role == sbol3.SBO_PRODUCT:
                sign = '+'
            else:
                raise ValueError (f'Unexpected role in {interaction.identity}: {tyto.SBO.get_term_by_uri(role)}')
            return f'{sign} {rate}*' + '*'.join(reactants)
        else:
            raise ValueError(f'No model for cleavage {interaction.name} in {interaction.identity}')
    elif i_type == sbol3.SBO_DEGRADATION:
        if len(interaction.participations) != 1:
            raise ValueError(f'Degradation assumed to have 1 participant, found {len(interaction.participations)} in {interaction.identity}')
        if f_type == sbol3.SBO_RNA:
            deg_rate = parameters[f'delta_g']
        elif f_type == sbol3.SBO_PROTEIN:
            species = variables[feature]
            deg_rate = parameters[f'delta_{species}']
        else:
            deg_rate = parameters[interaction]  # TODO: move this into actual parameters rather than name

        return f'- {deg_rate}*{variables[feature]}'

        pass
    elif i_type == sbol3.SBO_NON_COVALENT_BINDING:
        reactants = [variables[f] for f in all_in_role(interaction, sbol3.SBO_REACTANT)]
        [variables[f] for f in all_in_role(interaction, sbol3.SBO_PRODUCT)]  # Get products into the variable table
        rate = parameters[interaction]  # TODO: move this into actual parameters rather than name
        if role == sbol3.SBO_REACTANT:
            sign = '-'
        elif role == sbol3.SBO_PRODUCT:
            sign = '+'
        else:
            raise ValueError(f'Cannot handle type {tyto.SBO.get_term_by_uri(f_type)} in {interaction.identity}')
        return f'{sign} {rate}*' + '*'.join(reactants)
    elif i_type == sbol3.SBO_INHIBITION or i_type == sbol3.SBO_STIMULATION:
        # Pass for the regulation interactions that are taken care of in the other function, so you don't get a warning
        pass
    elif i_type == RECOMBINATION:
        if role == sbol3.SBO_MODIFIER or role == sbol3.SBO_MODIFIED:
            return None  # no effect on Cre concentration, not modeling excised element
        reactant = in_role(interaction, sbol3.SBO_REACTANT)
        recombinase = variables[in_role(interaction, sbol3.SBO_MODIFIER)]
        ct = containers[reactant]
        if len(ct) != 1:
            raise ValueError(f'Recombination expected 1 context, got {len(ct)} in {interaction.identity}')
        context = ct[0]
        rate = parameters['k_cre'] # TODO: move this into actual parameters rather than name
        if role == sbol3.SBO_REACTANT:
            sign = '-'
        elif role == sbol3.SBO_PRODUCT:
            sign = '+'
        else:
            raise ValueError(f'Cannot handle type {tyto.SBO.get_term_by_uri(f_type)} in {interaction.identity}')
        return f'{sign} {rate}*{variables[reactant]}*{recombinase}^4 + ' \
               f'({variables[feature]}/{variables[context]})*{differential(context)}'
    else:
        logging.warning(f'Cannot serialize interaction {interaction.identity} of type {tyto.SBO.get_term_by_uri(i_type)}')
        return None

# TODO: consider switch from ode45 to ode15s
ode_template = '''function [time_interval, y_out, y] = {}(time_span, parameters, initial, step)
% time_span is the hours values [start, stop]
% parameters is a Map of names to numbers (e.g., rate constants, decay rates, Hill coefficients)
% initial is a Map of variable names to initial values
% step is the number of hours between samples in output; defaults to 1
% Returns vector of time, matrix of output levels at those time points, matrix of all species
    if nargin < 4, step = 1; end
    
    % Define names for input/output variable indexes
    {}

    % Set initial values
    y0=zeros(1,{});
    {}
    
    % Run ODE
    solution = {}(@(t,x) diff_eq(t, x, parameters), time_span, y0);
    
    % Evaluate species levels at given times
    time_interval = time_span(1):step:time_span(end);
    y = deval(solution, time_interval);
    y_out = y([{}],:);
end

% ODE differential function
function dx=diff_eq(t, x, parameters)
    % Unpack parameters from parameter map
    {}
    
    % Unpack individual species from x
    x = max(1e-12,real(x)); % Truncate values just above zero
    {}
    
    % Compute derivative for each species
    {}
    
    % Pack derivatives for return, ensuring none are complex
    dx = real([{}])';
end
'''
"""Template for the Matlab simulation, including both the runner and the step function.
Format parameters are:

 1 protocol name
 2 Input/output variable names: VARIABLE = i
 3 Number of variables (integer)
 4 Initial value assignments for input variables: y0(VARIABLE) = initial(i)
 5 ODE function (ode45 or ode15s)
 6 Output indices: VARIABLE, VARIABLE, ...
 7 Parameter names: PARAMETER = i
 8 Unpacking of variables from x value: VARIABLE = x(i)
 9 Derivative equations for each species: dVARIABLE = EXPRESSION
 10 Packing of derivatives for return value: dVARIABLE, dVARIABLE, ...
"""


def format_model(name: str, parameters: List[str], variables: List[str], inputs: List[str], outputs: List[str],
                 derivatives: List[str], ode: str='ode45') -> str:
    """Generate a Matlab ODE simulation from the provided inputs

    :param name: protocol name
    :param parameters: list of parameter names
    :param variables: list of variable names
    :param inputs: list of names of input variables
    :param outputs: list of names of output variables
    :param derivatives: list of Matlab equations expressing the derivative for each variable
    :param ode: Matlab ODE function to use, defaults to ode45
    :return: string containing contents for Matlab simulation file
    """
    # Make the substructures
    parameter_names = "\n\t".join(f'{p} = parameters(\'{p}\');' for p in parameters)
    io_variable_names = "\n\t".join(f'{v} = {i};' for v, i in zip(variables, range(1, len(variables) + 1))
                                    if v in (set(inputs) | set(outputs)))
    initializations = "\n\t".join(f'y0({v}) = initial(\'{v}\');' for v in inputs)
    unpack_variables = "\n\t".join(f'{v} = x({i});' for v, i in zip(variables, range(1, len(variables) + 1)))
    pack_derivatives = ", ".join(f'{differential(v)}' for v in variables)
    return ode_template.format(name, io_variable_names, len(variables), initializations, ode, ", ".join(outputs),
                               parameter_names, unpack_variables, "\n\t".join(derivatives), pack_derivatives)


def make_matlab_model(system: sbol3.Component, ode: str='ode45') -> Tuple[str, List[str]]:
    """Generate a set of LaTeX equations for the identified system:

    :param system: system for which a model is to be generated
    :param ode: Matlab ODE function to use, defaults to ode45
    :return: string serialization of LaTeX equation collection
    """
    # for each feature, collect all of the interactions and constraints that it participates in
    interactions = {f: [i for i in system.interactions if [p for p in i.participations if p.participant == f.identity]]
                    for f in system.features}
    containers = {f: [c.subject.lookup() for c in system.constraints if c.restriction == sbol3.SBOL_CONTAINS and c.object == f.identity]
                  for f in system.features}
    regulators = {f: [c.subject.lookup() for c in system.constraints if c.restriction == sbol3.SBOL_MEETS and c.object == f.identity]
                  for f in system.features}
    regulation = {f: list(itertools.chain(*(interactions[r] for r in regulators[f]))) for f in regulators}

    # generate an ODE based on the roles in the interactions
    parameters = ParameterDictionary()  # dictionary of Interaction/string : parameter_name
    variables = VariableDictionary()  # dictionary of Feature: variable_name
    derivatives = []

    terms_added = set()
    for f in id_sort(system.features):
        interaction_terms = [t for t in [interaction_to_term(f, i, regulation, containers, parameters, variables)
                                         for i in id_sort(interactions[f])] if t]
        # If there is at least one term, then add an equation
        if interaction_terms:
            terms_added.add(f)
            derivatives.append(f'{differential(f)} = {" ".join(sorted(interaction_terms)).removeprefix("+")};')

    missing_terms = variables.keys() - terms_added
    derivatives = [f'{differential(f)} = 0;' for f in missing_terms] + derivatives

    # TODO: add d_VAR = 0 equations for any variables that didn't get an interaction term

    # Generate the actual document
    # TODO: interfaces will change to interface after resolution of https://github.com/SynBioDex/pySBOL3/issues/316
    # TODO: input/ouput will change to plural after resolution of https://github.com/SynBioDex/pySBOL3/issues/315
    parameter_names = sorted(set(parameters.values()))
    variable_names = sorted(variables.values())
    inputs = sorted([v for k, v in variables.items() if k.identity in (str(x) for x in system.interface.inputs)])
    outputs = sorted([v for k, v in variables.items() if k.identity in (str(x) for x in system.interface.outputs)])
    model = format_model(system.display_id, parameter_names, variable_names, inputs, outputs, derivatives, ode)

    return model, parameter_names
