import logging

import sbol3

from sbol_utilities.component import contains

#########################
# This file contains various patches and workarounds that will be deleted when their issues are resolved


def add_subfeature(container: sbol3.Feature, feature: sbol3.Feature) -> sbol3.Feature:
    """Deprecation wrapper for contained_by"""
    logging.warning(f'add_subfeature is deprecated; switch to contains')
    return contains(container, feature)

# May be resolved by https://github.com/SynBioDex/pySBOL3/issues/320
# If this is not the case, migrate up to SBOL-utilities
def identity_lt(self: sbol3.Identified, other: sbol3.Identified):
    return self.identity < other.identity
sbol3.Identified.__lt__ = identity_lt

