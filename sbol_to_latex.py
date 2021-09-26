import os

import sbol3

import latex_generation
from shared_global_names import *


doc = sbol3.Document()
doc.read(MODEL_FILE)

# For each system in the document, generate a matlab model
with open(os.path.join('generated_models', 'generated_equations.tex'), 'w') as out:
    for c in (o for o in doc.objects if isinstance(o, sbol3.Component)):
        out.write(latex_generation.make_latex_model(c))
