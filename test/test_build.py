import difflib
import unittest
import io

import sbol3

import builders
import latex_generation
from sbol_utilities.component import add_feature, constitutive, regulate, contains, add_interaction


class TestCircuitBuilding(unittest.TestCase):

    def test_kill_switch(self):
        """Make sure that the basic kill switch generates the right structure and from it the right LaTeX"""
        doc = sbol3.Document()
        sbol3.set_namespace('http://bbn.com/crispr-kill-switch/')

        system = sbol3.Component('Basic_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Basic Kill Switch")
        doc.add(system)

        aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
        sgRNA1_dna, genome = builders.make_crispr_module(aav)
        constitutive(sgRNA1_dna)
        # TODO: Warning will go away after resolution of hhttps://github.com/SynBioDex/pySBOL3/issues/324
        system.interface = sbol3.Interface(inputs=[aav, genome], outputs=[aav])

        generated = latex_generation.make_latex_model(system)
        expected = '''\\subsection{Basic Kill Switch}
\\label{s:Basic_kill_switch}
% Equations generated from http://bbn.com/crispr-kill-switch/Basic_kill_switch

\\begin{align}
\\diff{\\vectorGen{}}{t} & = - \\casCutRate{{}}\\vectorGen{}\\conc{\\cplx{\\cas}{1}}\\\\
\\diff{\\conc{\\cplx{\\cas}{1}}}{t} & =  \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{1}} - \\casCompDegradeRate{}\\conc{\\cplx{\\cas}{1}} - \\casCutRate{{}}\\vectorGen{}\\conc{\\cplx{\\cas}{1}}\\\\
\\diff{\\conc{\\cplx{\\cas}{2}}}{t} & =  \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{2}} - \\casCompDegradeRate{}\\conc{\\cplx{\\cas}{2}} - \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{}\\\\
\\diff{\\hostGen{}}{t} & = - \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{}\\\\
\\diff{\\conc{\\bound{\\cplx{\\cas}{1}}}}{t} & =  \\casCutRate{{}}\\vectorGen{}\\conc{\\cplx{\\cas}{1}} - \\casCompDegradeRate{}\\conc{\\bound{\\cplx{\\cas}{1}}}\\\\
\\diff{\\conc{\\bound{\\cplx{\\cas}{2}}}}{t} & =  \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{} - \\casCompDegradeRate{}\\conc{\\bound{\\cplx{\\cas}{2}}}\\\\
\\diff{\\edited{\\hostGen{}}}{t} & =  \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{}\\\\
\\diff{\\conc{\\proSp{\\cas{}}}}{t} & =  \\txtlRate{\\proSp{\\cas{}}}\\vectorGen{} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{1}} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{2}} - \\proDegradeRate{\\proSp{\\cas{}}}\\conc{\\proSp{\\cas{}}}\\\\
\\diff{\\conc{\\gRna{1}}}{t} & =  \\txRate{\\gRna{1}}\\vectorGen{} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{1}} - \\rnaDegradeRate{}\\conc{\\gRna{1}}\\\\
\\diff{\\conc{\\gRna{2}}}{t} & =  \\txRate{\\gRna{2}}\\vectorGen{} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{2}} - \\rnaDegradeRate{}\\conc{\\gRna{2}}
\\end{align}

'''

        diff = ''.join(difflib.unified_diff(io.StringIO(generated).readlines(), io.StringIO(expected).readlines(),
                                            fromfile='Generated', tofile='Expected'))
        assert not diff, f'Generated value does not match expectation: {diff}'

    def test_basic_tf_module(self):
        """Make sure that the basic TF module generates the right structure and from it the right LaTeX"""
        doc = sbol3.Document()
        sbol3.set_namespace('http://bbn.com/crispr-kill-switch/')
        system = sbol3.Component('Repression', sbol3.SBO_FUNCTIONAL_ENTITY, name="Simple Repression")
        doc.add(system)
        aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
        cas9_cds = contains(aav, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[sbol3.SO_CDS], name="Cas9-coding"))
        cas9 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name="Cas9"))
        add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {cas9_cds: sbol3.SBO_TEMPLATE, cas9: sbol3.SBO_PRODUCT})
        add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={cas9: sbol3.SBO_REACTANT})
        tf_cds, tf_promoter = builders.make_tf_module(aav, True)
        regulate(tf_promoter, cas9_cds)
        constitutive(tf_cds)

        generated = latex_generation.make_latex_model(system)
        expected = '''\\subsection{Simple Repression}
\\label{s:Repression}
% Equations generated from http://bbn.com/crispr-kill-switch/Repression

\\begin{align}
\\diff{\\conc{\\proSp{\\cas{}}}}{t} & =  \\txtlRate{\\proSp{\\cas{}}}\\frac{(K_R)^n}{(K_R)^n + \\conc{\\proSp{TF}}^n}\\vectorGen{} - \\proDegradeRate{\\proSp{\\cas{}}}\\conc{\\proSp{\\cas{}}}\\\\
\\diff{\\conc{\\proSp{TF}}}{t} & =  \\txtlRate{\\proSp{TF}}\\vectorGen{} - \\proDegradeRate{\\proSp{TF}}\\conc{\\proSp{TF}}
\\end{align}

'''

        diff = ''.join(difflib.unified_diff(io.StringIO(generated).readlines(), io.StringIO(expected).readlines(),
                                            fromfile='Generated', tofile='Expected'))
        assert not diff, f'Generated value does not match expectation: {diff}'

    def test_tf_kill_module(self):
        """Make sure that the TF-regulated kill switch generates the right structure and from it the right LaTeX"""
        doc = sbol3.Document()
        sbol3.set_namespace('http://bbn.com/crispr-kill-switch/')

        system = sbol3.Component('TF_delayed_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="TF Kill Switch")
        doc.add(system)
        aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
        sgRNA1_dna, genome = builders.make_crispr_module(aav)
        tf_cds, tf_promoter = builders.make_tf_module(aav, False)
        regulate(tf_promoter, sgRNA1_dna)
        constitutive(tf_cds)

        generated = latex_generation.make_latex_model(system)
        expected = '''\\subsection{TF Kill Switch}
\\label{s:TF_delayed_kill_switch}
% Equations generated from http://bbn.com/crispr-kill-switch/TF_delayed_kill_switch

\\begin{align}
\\diff{\\vectorGen{}}{t} & = - \\casCutRate{{}}\\vectorGen{}\\conc{\\cplx{\\cas}{1}}\\\\
\\diff{\\conc{\\cplx{\\cas}{1}}}{t} & =  \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{1}} - \\casCompDegradeRate{}\\conc{\\cplx{\\cas}{1}} - \\casCutRate{{}}\\vectorGen{}\\conc{\\cplx{\\cas}{1}}\\\\
\\diff{\\conc{\\cplx{\\cas}{2}}}{t} & =  \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{2}} - \\casCompDegradeRate{}\\conc{\\cplx{\\cas}{2}} - \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{}\\\\
\\diff{\\hostGen{}}{t} & = - \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{}\\\\
\\diff{\\conc{\\bound{\\cplx{\\cas}{1}}}}{t} & =  \\casCutRate{{}}\\vectorGen{}\\conc{\\cplx{\\cas}{1}} - \\casCompDegradeRate{}\\conc{\\bound{\\cplx{\\cas}{1}}}\\\\
\\diff{\\conc{\\bound{\\cplx{\\cas}{2}}}}{t} & =  \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{} - \\casCompDegradeRate{}\\conc{\\bound{\\cplx{\\cas}{2}}}\\\\
\\diff{\\edited{\\hostGen{}}}{t} & =  \\casCutRate{{}}\\conc{\\cplx{\\cas}{2}}\\hostGen{}\\\\
\\diff{\\conc{\\proSp{TF}}}{t} & =  \\txtlRate{\\proSp{TF}}\\vectorGen{} - \\proDegradeRate{\\proSp{TF}}\\conc{\\proSp{TF}}\\\\
\\diff{\\conc{\\proSp{\\cas{}}}}{t} & =  \\txtlRate{\\proSp{\\cas{}}}\\vectorGen{} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{1}} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{2}} - \\proDegradeRate{\\proSp{\\cas{}}}\\conc{\\proSp{\\cas{}}}\\\\
\\diff{\\conc{\\gRna{1}}}{t} & =  \\txRate{\\gRna{1}}\\frac{\\conc{\\proSp{TF}}^n}{(K_A)^n + \\conc{\\proSp{TF}}^n}\\vectorGen{} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{1}} - \\rnaDegradeRate{}\\conc{\\gRna{1}}\\\\
\\diff{\\conc{\\gRna{2}}}{t} & =  \\txRate{\\gRna{2}}\\vectorGen{} - \\gRnaBind{}\\conc{\\proSp{\\cas{}}}\\conc{\\gRna{2}} - \\rnaDegradeRate{}\\conc{\\gRna{2}}
\\end{align}

'''

        diff = ''.join(difflib.unified_diff(io.StringIO(generated).readlines(), io.StringIO(expected).readlines(),
                                            fromfile='Generated', tofile='Expected'))
        assert not diff, f'Generated value does not match expectation: {diff}'

    def test_basic_recombinase_module(self):
        """Make sure that the basic TF module generates the right structure and from it the right LaTeX"""
        doc = sbol3.Document()
        sbol3.set_namespace('http://bbn.com/crispr-kill-switch/')
        system = sbol3.Component('Recombinase', sbol3.SBO_FUNCTIONAL_ENTITY, name="Simple Recombinase")
        doc.add(system)
        aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
        cas9_cds = contains(aav, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[sbol3.SO_CDS], name="Cas9-coding"))
        cas9 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name="Cas9"))
        add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {cas9_cds: sbol3.SBO_TEMPLATE, cas9: sbol3.SBO_PRODUCT})
        add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={cas9: sbol3.SBO_REACTANT})
        cre_cds, cre_region = builders.make_recombinase_module(aav, True)
        regulate(cre_region, cas9_cds)
        constitutive(cre_cds)

        generated = latex_generation.make_latex_model(system)
        expected = '''\\subsection{Simple Recombinase}
\\label{s:Recombinase}
% Equations generated from http://bbn.com/crispr-kill-switch/Recombinase

\\begin{align}
\\diff{\\conc{\\proSp{\\cas{}}}}{t} & =  \\txtlRate{\\proSp{\\cas{}}}\\frac{\\edited{\\vectorGen{}_C}}{\\vectorGen{}}\\vectorGen{} - \\proDegradeRate{\\proSp{\\cas{}}}\\conc{\\proSp{\\cas{}}}\\\\
\\diff{\\conc{\\proSp{\cre{}}}}{t} & =  \\txtlRate{\\proSp{\\cre{}}}\\vectorGen{} - \\proDegradeRate{\\proSp{\\cre{}}}\\conc{\\proSp{\\cre{}}}\\\\
\\diff{\\vectorGen{}_C}{t} & = - \\creCutRate{} \\vectorGen{}_C \\conc{\\proSp{\\cre{}}}^4 + \\frac{\\vectorGen{}_C}{\\vectorGen{}} \diff{\\vectorGen{}}{t}\\\\
\\diff{\\edited{\\vectorGen{}_C}}{t} & =  \\creCutRate{} \\vectorGen{}_C \\conc{\\proSp{\\cre{}}}^4 + \\frac{\\edited{\\vectorGen{}_C}}{\\vectorGen{}} \\diff{\\vectorGen{}}{t}
\\end{align}

'''

        diff = ''.join(difflib.unified_diff(io.StringIO(generated).readlines(), io.StringIO(expected).readlines(),
                                            fromfile='Generated', tofile='Expected'))
        assert not diff, f'Generated value does not match expectation: {diff}'

    def test_two_tf_module(self):
        """Make sure that the basic TF module generates the right structure and from it the right LaTeX"""
        doc = sbol3.Document()
        sbol3.set_namespace('http://bbn.com/crispr-kill-switch/')
        system = sbol3.Component('TFTF', sbol3.SBO_FUNCTIONAL_ENTITY, name="Dual TF")
        doc.add(system)
        aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
        cas9_cds = contains(aav, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[sbol3.SO_CDS], name="Cas9-coding"))
        cas9 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name="Cas9"))
        add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {cas9_cds: sbol3.SBO_TEMPLATE, cas9: sbol3.SBO_PRODUCT})
        add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={cas9: sbol3.SBO_REACTANT})
        tf_cds, tf_target = builders.make_tf_module(aav, True)
        tf2_cds, tf2_target = builders.make_tf_module(aav, False, True)
        regulate(tf_target, cas9_cds)
        regulate(tf2_target, tf_cds)
        constitutive(tf2_cds)

        generated = latex_generation.make_latex_model(system)
        expected = '''\\subsection{Dual TF}
\\label{s:TFTF}
% Equations generated from http://bbn.com/crispr-kill-switch/TFTF

\\begin{align}
\\diff{\\conc{\\proSp{\\cas{}}}}{t} & =  \\txtlRate{\\proSp{\\cas{}}}\\frac{(K_R)^n}{(K_R)^n + \\conc{\\proSp{TF}}^n}\\vectorGen{} - \\proDegradeRate{\\proSp{\\cas{}}}\\conc{\\proSp{\\cas{}}}\\\\
\\diff{\\conc{\\proSp{TF}}}{t} & =  \\txtlRate{\\proSp{TF}}\\frac{\\conc{\\proSp{TF2}}^n}{(K_A)^n + \\conc{\\proSp{TF2}}^n}\\vectorGen{} - \\proDegradeRate{\\proSp{TF}}\\conc{\\proSp{TF}}\\\\
\\diff{\\conc{\\proSp{TF2}}}{t} & =  \\txtlRate{\\proSp{TF2}}\\vectorGen{} - \\proDegradeRate{\\proSp{TF2}}\\conc{\\proSp{TF2}}
\\end{align}

'''

        diff = ''.join(difflib.unified_diff(io.StringIO(generated).readlines(), io.StringIO(expected).readlines(),
                                            fromfile='Generated', tofile='Expected'))
        assert not diff, f'Generated value does not match expectation: {diff}'

    def test_two_recombinase_module(self):
        """Make sure that the basic TF module generates the right structure and from it the right LaTeX"""
        doc = sbol3.Document()
        sbol3.set_namespace('http://bbn.com/crispr-kill-switch/')
        system = sbol3.Component('Recombinase', sbol3.SBO_FUNCTIONAL_ENTITY, name="Dual Recombinase")
        doc.add(system)
        aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
        cas9_cds = contains(aav, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[sbol3.SO_CDS], name="Cas9-coding"))
        cas9 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name="Cas9"))
        add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {cas9_cds: sbol3.SBO_TEMPLATE, cas9: sbol3.SBO_PRODUCT})
        add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={cas9: sbol3.SBO_REACTANT})
        cre_cds, cre_region = builders.make_recombinase_module(aav, True)
        creH_cds, creH_region = builders.make_recombinase_module(aav, False, True)
        regulate(cre_region, cas9_cds)
        regulate(creH_region, cre_cds)
        constitutive(creH_cds)

        generated = latex_generation.make_latex_model(system)
        expected = '''\\subsection{Dual Recombinase}
\\label{s:Recombinase}
% Equations generated from http://bbn.com/crispr-kill-switch/Recombinase

\\begin{align}
\\diff{\\conc{\\proSp{CreH}}}{t} & =  \\txtlRate{\\proSp{CreH}}\\vectorGen{} - \\proDegradeRate{\\proSp{CreH}}\\conc{\\proSp{CreH}}\\\\
\\diff{\\vectorGen{}_{CH}}{t} & = - \\creCutRate{} \\vectorGen{}_{CH} \\conc{\\proSp{CreH}}^4 + \\frac{\\vectorGen{}_{CH}}{\\vectorGen{}} \\diff{\\vectorGen{}}{t}\\\\
\\diff{\\edited{\\vectorGen{}_{CH}}}{t} & =  \\creCutRate{} \\vectorGen{}_{CH} \\conc{\\proSp{CreH}}^4 + \\frac{\\edited{\\vectorGen{}_{CH}}}{\\vectorGen{}} \\diff{\\vectorGen{}}{t}\\\\
\\diff{\\conc{\\proSp{\\cas{}}}}{t} & =  \\txtlRate{\\proSp{\\cas{}}}\\frac{\\edited{\\vectorGen{}_C}}{\\vectorGen{}}\\vectorGen{} - \\proDegradeRate{\\proSp{\\cas{}}}\\conc{\\proSp{\\cas{}}}\\\\
\\diff{\\conc{\\proSp{\\cre{}}}}{t} & =  \\txtlRate{\\proSp{\\cre{}}}\\frac{\\vectorGen{}_{CH}}{\\vectorGen{}}\\vectorGen{} - \\proDegradeRate{\\proSp{\\cre{}}}\\conc{\\proSp{\\cre{}}}\\\\
\\diff{\\vectorGen{}_C}{t} & = - \\creCutRate{} \\vectorGen{}_C \\conc{\\proSp{\\cre{}}}^4 + \\frac{\\vectorGen{}_C}{\\vectorGen{}} \diff{\\vectorGen{}}{t}\\\\
\\diff{\\edited{\\vectorGen{}_C}}{t} & =  \\creCutRate{} \\vectorGen{}_C \\conc{\\proSp{\\cre{}}}^4 + \\frac{\\edited{\\vectorGen{}_C}}{\\vectorGen{}} \\diff{\\vectorGen{}}{t}
\\end{align}

'''

        diff = ''.join(difflib.unified_diff(io.StringIO(generated).readlines(), io.StringIO(expected).readlines(),
                                            fromfile='Generated', tofile='Expected'))
        assert not diff, f'Generated value does not match expectation: {diff}'

if __name__ == '__main__':
    unittest.main()
