#!/usr/bin/python
import nest
import nest.voltage_trace
import pylab as pl

# Install our own module
nest.Install('coronetmodule')

# always reset the Kernel
nest.ResetKernel()

Neuron = nest.Create('coronet_neuron')
PoissonStim = nest.Create('poisson_generator', 2, params={'rate':20.})

nest.Connect([PoissonStim[0]], Neuron, params={'weight':20.})
nest.Connect([PoissonStim[1]], Neuron, params={'weight':-10.})

vm = nest.Create('voltmeter')
nest.Connect(vm,Neuron)

# recording conductances:
mm = nest.Create('multimeter', params={'record_from':['g_ex','g_in'], 'interval':0.1})
nest.Connect(mm,Neuron)


nest.Simulate(1000.)

pl.subplot(211)
nest.voltage_trace.from_device(vm)

pl.subplot(212)
recorded_events = nest.GetStatus(mm)[0]['events']
pl.plot(recorded_events['times'], recorded_events['g_ex'], label='g_ex')
pl.plot(recorded_events['times'], recorded_events['g_in'], label='g_in')
pl.legend()
pl.xlabel("Time [ms]")
pl.ylabel("Synaptic Conductance [nS]")
pl.show()
