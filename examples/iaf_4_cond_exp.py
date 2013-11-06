#!/usr/bin/python
import nest
import nest.voltage_trace
import pylab as pl

# Install our own module
nest.Install('coronetmodule')

# always reset the Kernel
nest.ResetKernel()

Neuron = nest.Create('iaf_4_cond_exp')
PoissonStim = nest.Create('poisson_generator', 2, params={'rate':20.})

# Get dictionary of receptor types for coronet neuron
receptors = nest.GetDefaults('iaf_4_cond_exp')['receptor_types']

# before (fails now!)
#nest.Connect([PoissonStim[0]], Neuron, params={'weight':20.})
#nest.Connect([PoissonStim[1]], Neuron, params={'weight':-10.})

# now with receptor_type:
nest.Connect([PoissonStim[0]], Neuron, params={'weight':20., 'receptor_type':receptors['SYN_1']})
nest.Connect([PoissonStim[1]], Neuron, params={'weight':10., 'receptor_type':receptors['SYN_2']})

vm = nest.Create('voltmeter')
nest.Connect(vm,Neuron)

# recording conductances:
mm = nest.Create('multimeter', params={'record_from':['g_syn_1','g_syn_2'], 'interval':0.1})
nest.Connect(mm,Neuron)


nest.Simulate(1000.)
print "Done"

pl.subplot(211)
nest.voltage_trace.from_device(vm)

pl.subplot(212)
recorded_events = nest.GetStatus(mm)[0]['events']
pl.plot(recorded_events['times'], recorded_events['g_syn_1'], label='g_syn_1')
pl.plot(recorded_events['times'], recorded_events['g_syn_2'], label='g_syn_2')
pl.legend()
pl.xlabel("Time [ms]")
pl.ylabel("Synaptic Conductance [nS]")
pl.show()
