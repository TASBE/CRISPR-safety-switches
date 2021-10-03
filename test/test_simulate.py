import filecmp
import math
import os
import tempfile
import unittest
from shutil import copy

import numpy
import sbol3
import oct2py

import builders
import matlab_generation
from sbol_utilities.component import add_feature, contains, add_interaction, regulate, constitutive


class TestSimulation(unittest.TestCase):

    def test_kill_switch(self):
        """Make sure running the basic kill switch produces reasonable values"""

        # build the kill switch
        doc = sbol3.Document()
        sbol3.set_namespace('http://bbn.com/crispr-kill-switch/')

        system = sbol3.Component('Basic_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Basic Kill Switch")
        doc.add(system)
        aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
        sgRNA1_dna, genome = builders.make_crispr_module(aav)
        builders.constitutive(sgRNA1_dna)
        # TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/315
        # TODO: interfaces will change to interface after resolution of https://github.com/SynBioDex/pySBOL3/issues/316
        system.interfaces = sbol3.Interface(input=[aav, genome], output=[aav])

        # generate the matlab and write it to a temp file
        tmp_dir = tempfile.mkdtemp()
        test_dir = os.path.dirname(os.path.realpath(__file__))
        copy(os.path.join(test_dir, 'test_files', 'deval_octave.m'), os.path.join(tmp_dir, 'deval.m'))
        with open(os.path.join(tmp_dir, 'Basic_kill_switch.m'), 'w') as f:
            model, parameters = matlab_generation.make_matlab_model(system)
            f.write(model)

        # check the the model generated is as expected
        comparison_file = os.path.join(test_dir, 'test_files', 'Basic_kill_switch.m')
        assert filecmp.cmp(os.path.join(tmp_dir, 'Basic_kill_switch.m'), comparison_file)

        # run the simulation and make sure outcomes are reasonable
        oc = oct2py.Oct2Py()
        oc.eval(f'cd {tmp_dir}')
        oc.eval(f'parameters = containers.Map();')
        for p in parameters:
            oc.eval(f'parameters(\'{p}\') = 1;')
        oc.eval('initial = containers.Map();')
        oc.eval('initial(\'AAV\') = 10;')
        oc.eval('initial(\'genome\') = 1;')

        oc.eval('[t,y_out,y] = Basic_kill_switch([0 72], parameters, initial);')
        # check variable values in numpy arrays
        t = oc.pull('t')
        assert isinstance(t, numpy.ndarray)
        # time progresses reasonably
        assert t.size == 73 and t[0, 0] == 0 and t[0, -1] == 72
        # 10 plasmids decay to 1 over time
        y_out = oc.pull('y_out')
        assert isinstance(y_out, numpy.ndarray)
        assert y_out.size == t.size
        assert y_out[0, 0] > 9.99999
        assert y_out[0, -1] < 0.15

    def test_basic_tf_module(self):
        """Make sure that the basic TF module generates the right structure and from it the right LaTeX"""
        doc = sbol3.Document()
        sbol3.set_namespace('http://bbn.com/crispr-kill-switch/')
        system = sbol3.Component('simple_repression', sbol3.SBO_FUNCTIONAL_ENTITY, name="Simple Repression")
        doc.add(system)
        aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
        cas9_cds = contains(aav, sbol3.LocalSubComponent([sbol3.SBO_DNA], roles=[sbol3.SO_CDS], name="Cas9-coding"))
        cas9 = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_PROTEIN], name="Cas9"))
        add_interaction(sbol3.SBO_GENETIC_PRODUCTION, {cas9_cds: sbol3.SBO_TEMPLATE, cas9: sbol3.SBO_PRODUCT})
        add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={cas9: sbol3.SBO_REACTANT})
        tf_cds, tf_promoter = builders.make_tf_module(aav, True)
        regulate(tf_promoter, cas9_cds)
        constitutive(tf_cds)
        # TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/315
        # TODO: interfaces will change to interface after resolution of https://github.com/SynBioDex/pySBOL3/issues/316
        system.interfaces = sbol3.Interface(input=[aav], output=[cas9])

        # generate the matlab and write it to a temp file
        tmp_dir = tempfile.mkdtemp()
        test_dir = os.path.dirname(os.path.realpath(__file__))
        copy(os.path.join(test_dir, 'test_files', 'deval_octave.m'), os.path.join(tmp_dir, 'deval.m'))
        with open(os.path.join(tmp_dir, 'simple_repression.m'), 'w') as f:
            model, parameters = matlab_generation.make_matlab_model(system)
            f.write(model)

        # check the the model generated is as expected
        comparison_file = os.path.join(test_dir, 'test_files', 'simple_repression.m')
        assert filecmp.cmp(os.path.join(tmp_dir, 'simple_repression.m'), comparison_file)

        copy(os.path.join(test_dir, 'test_files', 'run_tf_sim.m'), tmp_dir)
        oc = oct2py.Oct2Py()
        oc.eval(f'cd {tmp_dir}')
        oc.eval('run_tf_sim')
        # check variable values in numpy arrays
        assert oc.eval('abs(test_end_points - expected) < 0.01;').all()

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
        add_interaction(sbol3.SBO_DEGRADATION, name='Cas degradation', participants={cas9: sbol3.SBO_REACTANT})
        cre_cds, cre_region = builders.make_recombinase_module(aav, True)
        regulate(cre_region, cas9_cds)
        constitutive(cre_cds)
        # TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/315
        # TODO: interfaces will change to interface after resolution of https://github.com/SynBioDex/pySBOL3/issues/316
        system.interfaces = sbol3.Interface(input=[aav, cre_region], output=[cas9])

        # generate the matlab and write it to a temp file
        tmp_dir = tempfile.mkdtemp()
        test_dir = os.path.dirname(os.path.realpath(__file__))
        copy(os.path.join(test_dir, 'test_files', 'deval_octave.m'), os.path.join(tmp_dir, 'deval.m'))
        with open(os.path.join(tmp_dir, 'simple_recombinase.m'), 'w') as f:
            model, parameters = matlab_generation.make_matlab_model(system)
            f.write(model)

        # check the the model generated is as expected
        comparison_file = os.path.join(test_dir, 'test_files', 'simple_recombinase.m')
        assert filecmp.cmp(os.path.join(tmp_dir, 'simple_recombinase.m'), comparison_file)

        copy(os.path.join(test_dir, 'test_files', 'run_cre_sim.m'), tmp_dir)
        oc = oct2py.Oct2Py()
        oc.eval(f'cd {tmp_dir}')
        oc.eval('run_cre_sim')
        # check variable values in numpy arrays
        assert oc.eval('abs(test_end_points - expected) < 0.01;').all()



if __name__ == '__main__':
    unittest.main()
