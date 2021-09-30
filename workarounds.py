import logging

import sbol3

from sbol_utilities.component import contains

#########################
# This file contains various patches and workarounds that will be deleted when their issues are resolved

def add_subfeature(container: sbol3.Feature, feature: sbol3.Feature) -> sbol3.Feature:
    "Deprecation wrapper for contained_by"
    logging.warning(f'add_subfeature is deprecated; switch to contains')
    return contains(container, feature)
