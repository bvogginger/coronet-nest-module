import unittest
import numpy as np

import nest
# Install our own module
nest.Install('coronetmodule')



class CoronetNeuronTest(unittest.TestCase):
    """
    Test for the NEST 'coronet_neuron'
    """
    def runTest(self):
        """
        we compare if the coronet neuron does the same like the iaf_cond_exp neuron
        """
        # always reset the Kernel
        nest.ResetKernel()

        cn = nest.Create('coronet_neuron')
        iaf = nest.Create('iaf_cond_exp')

        # create two stimuli via parrot neurons, such that both neurons receive the same stimulus
        PoissonStim = nest.Create('poisson_generator', 2, params={'rate':20.})
        ParrotStim = nest.Create('parrot_neuron', 2)
        nest.Connect(PoissonStim, ParrotStim)

        # Get dictionary of receptor types for coronet neuron
        receptors = nest.GetDefaults('coronet_neuron')['receptor_types']

        w_ex = 20.
        w_in = 10.

        # connect to iaf
        nest.Connect([ParrotStim[0]], iaf, params={'weight':w_ex})
        nest.Connect([ParrotStim[1]], iaf, params={'weight':-w_in})

        # now with receptor_type:
        nest.Connect([ParrotStim[0]], cn, params={'weight':w_ex, 'receptor_type':receptors['EX']})
        nest.Connect([ParrotStim[1]], cn, params={'weight':w_in, 'receptor_type':receptors['IN']})

        vm_cn = nest.Create('voltmeter')
        nest.Connect(vm_cn,cn)

        vm_iaf = nest.Create('voltmeter')
        nest.Connect(vm_iaf,iaf)

        nest.Simulate(1000.)

        vms_cn = nest.GetStatus(vm_cn)[0]['events']['V_m']
        vms_iaf = nest.GetStatus(vm_iaf)[0]['events']['V_m']
        self.assertTrue( np.all(vms_cn == vms_iaf), "Voltage Traces of coronet_neuron and iaf_cond_exp are not equal")

class IAF4CondExpNeuronTest(unittest.TestCase):
    """
    Test for the NEST 'iaf_4_cond_exp'
    """
    def runTest(self):
        """
        we compare if the iaf_4_cond_exp neuron does the same like the iaf_cond_exp neuron
        """
        # always reset the Kernel
        nest.ResetKernel()

        cn = nest.Create('iaf_4_cond_exp')
        iaf = nest.Create('iaf_cond_exp')

        # set reversal potentials
        iaf_defaults = nest.GetDefaults('iaf_cond_exp')
        nest.SetStatus(cn, {
            "E_syn_1":iaf_defaults["E_ex"],
            "tau_syn_1":iaf_defaults["tau_syn_ex"],
            "E_syn_2":iaf_defaults["E_in"],
            "tau_syn_2":iaf_defaults["tau_syn_in"]
            })

        # create two stimuli via parrot neurons, such that both neurons receive the same stimulus
        PoissonStim = nest.Create('poisson_generator', 2, params={'rate':20.})
        ParrotStim = nest.Create('parrot_neuron', 2)
        nest.Connect(PoissonStim, ParrotStim)

        # Get dictionary of receptor types for coronet neuron
        receptors = nest.GetDefaults('iaf_4_cond_exp')['receptor_types']
        print receptors

        w_ex = 20.
        w_in = 10.

        # connect to iaf
        nest.Connect([ParrotStim[0]], iaf, params={'weight':w_ex})
        nest.Connect([ParrotStim[1]], iaf, params={'weight':-w_in})

        # now with receptor_type:
        nest.Connect([ParrotStim[0]], cn, params={'weight':w_ex, 'receptor_type':receptors['SYN_1']})
        nest.Connect([ParrotStim[1]], cn, params={'weight':w_in, 'receptor_type':receptors['SYN_2']})

        vm_cn = nest.Create('voltmeter')
        nest.Connect(vm_cn,cn)

        vm_iaf = nest.Create('voltmeter')
        nest.Connect(vm_iaf,iaf)

        nest.Simulate(1000.)

        vms_cn = nest.GetStatus(vm_cn)[0]['events']['V_m']
        vms_iaf = nest.GetStatus(vm_iaf)[0]['events']['V_m']
        self.assertTrue( np.all(vms_cn == vms_iaf), "Voltage Traces of iaf_4_cond_exp and iaf_cond_exp are not equal")

if __name__ == "__main__":
    CoronetNeuronTest().runTest()
    IAF4CondExpNeuronTest().runTest()
