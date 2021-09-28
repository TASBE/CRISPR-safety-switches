from typing import Optional, Dict, List

import sbol3

###########################
# Issue workarounds


# TODO: remove after resolution of https://github.com/SynBioDex/pySBOL3/issues/234
from sbol_utilities.helper_functions import id_sort


def get_parent(self: sbol3.Identified) -> Optional[sbol3.Identified]:
    """Find the parent of this child object

    :param self: Object to search from
    :return: parent or None if it cannot be found (e.g., this is a toplevel)
    """
    if self.identity:
        return self.document.find(self.identity.rsplit('/', 1)[0])
    else:
        return None
sbol3.Identified.get_parent = get_parent


# TODO: remove after resolution of https://github.com/SynBioDex/pySBOL3/issues/234
def get_toplevel(self: sbol3.Identified) -> Optional[sbol3.TopLevel]:
    """Find the SBOL3 TopLevel object containing this SBOL3 object

    :param self: Object to search from
    :return: Enclosing TopLevel (self if it is a TopLevel) or None if there is nothing enclosing
    """
    if isinstance(self, sbol3.TopLevel):
        return self
    else:
        parent = self.get_parent()
        if parent:
            return get_toplevel(parent)
        else:
            return None


###########################
# Helper functions

# TODO: migrate these to sbol_utilities

def add_feature(system: sbol3.Component, feature: sbol3.Feature) -> sbol3.Feature:
    """Pass-through adder for allowing slightly more compact code

    :param system: system to add feature to
    :param feature: to be added to system
    :return: feature, returned without modification
    """
    system.features.append(feature)
    return feature


def add_subfeature(container: sbol3.Feature, feature: sbol3.Feature) -> sbol3.Feature:
    """Add a sequence feature within a larger containing sequence feature, along with the topological constraint

    :param container: feature in the system that will contain this feature
    :param feature: to be added to system
    :return: feature, returned without modification
    """
    system = get_parent(container)
    system.features.append(feature)
    c = sbol3.Constraint(sbol3.SBOL_CONTAINS, subject=container, object=feature)
    system.constraints.append(c)
    return feature


def add_interaction(system: sbol3.Component, interaction_type: str, participants: Dict[sbol3.Feature, str],
                    name: str = None) -> sbol3.Interaction:
    """Compact function for creation of an interaction

    :param system: system to add interaction to
    :param interaction_type: SBO type of interaction to be to be added
    :param participants: dictionary assigning features to roles for participations
    :param name: name for the interaction
    :return: interaction
    """
    participations = [sbol3.Participation([r], p) for p, r in participants.items()]
    interaction = sbol3.Interaction([interaction_type], participations=participations, name=name)
    system.interactions.append(interaction)
    return interaction


def regulate(five_prime: sbol3.Feature, target: sbol3.Feature) -> sbol3.Constraint:
    """Connect a 5' regulatory region to control the expression of a CDS or ncRNA target

    :param five_prime: Regulatory region to place upstream of target
    :param target: CDS or ncRNA region to be regulated
    :return: Meets constraint used to link the two
    """
    system = get_parent(five_prime)
    c = sbol3.Constraint(sbol3.SBOL_MEETS, subject=five_prime, object=target)
    system.constraints.append(c)
    return c


def in_role(interaction: sbol3.Interaction, role: str) -> sbol3.Feature:
    """Find the (precisely one) feature with a given role in the interaction

    :param interaction: interaction to search
    :param role: role to search for
    :return Feature playing that role
    """
    feature_participation = [p for p in interaction.participations if role in p.roles]
    if len(feature_participation) != 1:
        raise ValueError(f'Role can be in 1 participant: found {len(feature_participation)} in {interaction.identity}')
    return feature_participation[0].participant.lookup()


def all_in_role(interaction: sbol3.Interaction, role: str) -> List[sbol3.Feature]:
    """Find the features with a given role in the interaction

    :param interaction: interaction to search
    :param role: role to search for
    :return sorted list of Features playing that role
    """
    return id_sort([p.participant.lookup() for p in interaction.participations if role in p.roles])