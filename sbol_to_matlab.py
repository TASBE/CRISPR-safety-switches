import os

import sbol3

def make_matlab_model(system: sbol3.Component) -> str:
    """Generate a Matlab model for the identified system:

    :param system: system for which a model is to be generated
    :return: string serialization of Matlab model
    """
    # for each feature, collect all of the interactions that it participates in
    # generate an ODE based on the roles in the interactions

    # core simulation is just a list of all of the ODE lines, packed and unpacked

    # wrapper simulation


doc = sbol3.Document()
# For each system in the document, generate a matlab model

for c in (o for o in doc.objects if isinstance(o, sbol3.Component)):
    with open(os.path.join('generated_models', f'{c.display_id}.m'), 'w') as f:
        f.write(make_matlab_model(c))
