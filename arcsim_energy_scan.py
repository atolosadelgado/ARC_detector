# energy_list = [x for x in range(1,10,2)] + [x for x in range(10,60,20)]
energy_list = [2,5,7,10,15,50]
import os

multiplicity = 100
particle_list = [ 'e+', 'pi+', 'mu+', 'kaon+','proton']

for particle in particle_list:
   for energy_GeV in energy_list:
       outputFile = f'arcsim_{particle}_{energy_GeV}GeV.root'
       cmd = f"nohup ddsim --steeringFile steering.py --gun.energy '{energy_GeV}*GeV' --gun.multiplicity {multiplicity} --gun.particle '{particle}' --outputFile {outputFile} &"

       os.system( cmd )
