import filecmp
import os
import tempfile
import unittest
from shutil import copy

import numpy
import sbol3
import oct2py

import builders
import matlab_generation
from sbol_utilities.component import add_feature


class TestSimulation(unittest.TestCase):

    def test_kill_switch(self):
        """Make sure running the basic kill switch produces reasonable values"""

        # build the kill switch
        doc = sbol3.Document()
        sbol3.set_namespace('http://bbn.com/crispr-kill-switch/')

        system = sbol3.Component('Basic_kill_switch', sbol3.SBO_FUNCTIONAL_ENTITY, name="Basic Kill Switch")
        doc.add(system)
        aav = add_feature(system, sbol3.LocalSubComponent([sbol3.SBO_DNA], name='AAV'))
        sgRNA1_dna, genome, sgRNA1_rna = builders.make_crispr_module(aav)
        builders.constitutive(sgRNA1_dna)
        # TODO: Warning will go away after resolution of https://github.com/SynBioDex/pySBOL3/issues/315
        interface = sbol3.Interface(input=[aav, genome], output=[aav])
        # TODO: interfaces will change to interface after resolution of https://github.com/SynBioDex/pySBOL3/issues/316
        system.interfaces = interface

        # generate the matlab and write it to a temp file
        tmp_dir = tempfile.mkdtemp()
        test_dir = os.path.dirname(os.path.realpath(__file__))
        copy(os.path.join(test_dir, 'test_files', 'deval.m'), tmp_dir)
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
        oc.eval('[t,y_out,y] = Basic_kill_switch([0 72], parameters, [10 1]);')
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

if __name__ == '__main__':
    unittest.main()
