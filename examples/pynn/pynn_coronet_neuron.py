import pyNN.nest as sim
sim.nest.Install('coronetmodule')
sim.setup()
sim.nest.SetKernelStatus({'dict_miss_is_error': False})

coronet_neuron = sim.native_cell_type('coronet_neuron')
p1 = sim.Population(10, coronet_neuron)
s1 = sim.Population(10, sim.SpikeSourcePoisson, {'rate':20})
s2 = sim.Population(10, sim.SpikeSourcePoisson, {'rate':20})
sim.Projection(s1, p1, sim.OneToOneConnector(weights=0.02), target="EX")
sim.Projection(s2, p1, sim.OneToOneConnector(weights=0.01), target="IN")

p1.initialize('V_m', -76.) # this works as expected

p1.record()
p1._record('V_m') # ugly
sim.run(1000)
id, t, v = p1.recorders['V_m'].get().T # ugly

import pylab as pl
import numpy as np
pl.figure()
sl = p1.getSpikes()
pl.plot(sl[:,1],sl[:,0],'.')
pl.figure()
id_is_0 = np.where(id==0)
pl.plot(t[id_is_0],v[id_is_0])
pl.show()
