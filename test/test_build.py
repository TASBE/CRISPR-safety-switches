import difflib
import unittest
import io

import sbol3

import builders
import latex_generation
from helpers import add_feature


class test_build(unittest.TestCase):

    def test_kill_switch(self):
        """Make sure that the basic kill switch generates the right structure and from it the right LaTeX"""
        doc = sbol3.Document()
        sbol3.set_namespace('http://bbn.com/crispr-kill-switch/')

        system = sbol3.Component('Basic_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Basic Kill Switch")
        doc.add(system)

        aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
        sgRNA1_dna, genome = builders.make_crispr_module(aav)
        builders.constitutive(sgRNA1_dna)
        # TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/315
        interface = sbol3.Interface(input=[aav, genome], output=[aav])
        # TODO: interfaces will change to interface after resolution of https://github.com/SynBioDex/pySBOL3/issues/316
        system.interfaces = interface

        generated = latex_generation.make_latex_model(system)
        expected = '''\\subsection{Basic Kill Switch}
\\label{s:Basic_kill_switch}
% Equations generated from http://bbn.com/crispr-kill-switch/Basic_kill_switch

\\begin{align}
\\diff{\\vectorGen{}}{t} & = - \\casCutRate{{}}\\vectorGen{}\\conc{\\cplx{\\cas}{1}}\\\\ 
\\diff{\\conc{\\cplx{\\cas}{1}}}{t} & =  \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{1}} - \\casCutRate{{}}\\vectorGen{}\\conc{\\cplx{\\cas}{1}} - \\casCompDegradeRate{}\\conc{\\cplx{\\cas}{1}}\\\\ 
\\diff{\\conc{\\cplx{\\cas}{2}}}{t} & =  \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{2}} - \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{} - \\casCompDegradeRate{}\\conc{\\cplx{\\cas}{2}}\\\\ 
\\diff{\\hostGen{}}{t} & = - \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{}\\\\ 
\\diff{\\conc{\\bound{\\cplx{\\cas}{1}}}}{t} & = - \\casCompDegradeRate{}\\conc{\\bound{\\cplx{\\cas}{1}}} + \\casCutRate{{}}\\vectorGen{}\\conc{\\cplx{\\cas}{1}}\\\\ 
\\diff{\\conc{\\bound{\\cplx{\\cas}{2}}}}{t} & = - \\casCompDegradeRate{}\\conc{\\bound{\\cplx{\\cas}{2}}} + \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{}\\\\ 
\\diff{\\edited{\\hostGen{}}}{t} & =  \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{}\\\\ 
\\diff{\\conc{\\proSp{\\cas{}}}}{t} & =  \\txtlRate{\\proSp{\\cas{}}}\\vectorGen{} - \\proDegradeRate{\\proSp{\\cas{}}}\\conc{\\proSp{\\cas{}}} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{1}} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{2}}\\\\ 
\\diff{\\conc{\\gRna{1}}}{t} & =  \\txRate{\\gRna{1}}\\vectorGen{} - \\rnaDegradeRate{}\\conc{\\gRna{1}} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{1}}\\\\ 
\\diff{\\conc{\\gRna{2}}}{t} & =  \\txRate{\\gRna{2}}\\vectorGen{} - \\rnaDegradeRate{}\\conc{\\gRna{2}} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{2}}
\\end{align}

'''

        diff = ''.join(difflib.unified_diff(io.StringIO(generated).readlines(), io.StringIO(expected).readlines(),
                                            fromfile='Generated', tofile='Expected'))
        assert not diff, f'Generated value does not match expectation: {diff}'


if __name__ == '__main__':
    unittest.main()
