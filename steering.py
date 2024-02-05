#steering.py 
from DDSim.DD4hepSimulation import DD4hepSimulation
SIM = DD4hepSimulation()
from g4units import mm, GeV, MeV

# Register optical physics and Cerenkov process to the default physics
def setupCerenkov(kernel):
        from DDG4 import PhysicsList

        seq = kernel.physicsList()
        cerenkov = PhysicsList(kernel, "Geant4CerenkovPhysics/CerenkovPhys")
        cerenkov.MaxNumPhotonsPerStep = 10
        cerenkov.MaxBetaChangePerStep = 10.0
        cerenkov.TrackSecondariesFirst = False
        cerenkov.VerboseLevel = 0
        cerenkov.enableUI()
        seq.adopt(cerenkov)
        ph = PhysicsList(kernel, "Geant4OpticalPhotonPhysics/OpticalGammaPhys")
        ph.addParticleConstructor("G4OpticalPhoton")
        ph.VerboseLevel = 0
        ph.BoundaryInvokeSD = True
        ph.enableUI()
        seq.adopt(ph)
        return None
SIM.physics.setupUserPhysics(setupCerenkov)

# Associate the Geant4OpticalTrackerAction to these detectors
# this action register total energy of the opt photon as a single hit
# and kills the optical photon, so no time is wasted tracking them
SIM.action.mapActions["ARCBARREL"] = "Geant4OpticalTrackerAction"
SIM.action.mapActions["ARCENDCAP"] = "Geant4OpticalTrackerAction"
SIM.action.mapActions["ARC_DETECTORNAME"] = "Geant4OpticalTrackerAction"

# Disable user tracker particle handler, so hits can be associated to photons
SIM.part.userParticleHandler = ""

# Register hit with low energy, compatible with zero
SIM.filter.tracker = "edep0"

# Define filter, so detector is only sensitive to optical photons
SIM.filter.filters["opticalphotons"] = dict(
        name="ParticleSelectFilter/OpticalPhotonSelector",
        parameter={"particle": "opticalphoton"},
        )
SIM.filter.mapDetFilter["ARCBARREL"] = "opticalphotons"
SIM.filter.mapDetFilter["ARCENDCAP"] = "opticalphotons"
SIM.filter.mapDetFilter["ARC_DETECTORNAME"] = "opticalphotons"

# Particle gun settings: pions with fixed energy, random direction
SIM.numberOfEvents = 1
SIM.enableGun = True
SIM.gun.energy = "50*GeV"
SIM.gun.particle = "e+"
SIM.gun.multiplicity = 1
SIM.gun.position = "0 0 -20*cm"
SIM.gun.direction = "0 0 1"
# SIM.gun.momentumMin = 0.1*GeV
# SIM.gun.momentumMax = 0.11*GeV
SIM.compactFile = "./compact/arc_v0.xml"
SIM.outputFile = "arcsim_e+_50GeV_key4hep.root"
SIM.runType = "batch"
#SIM.runType = "qt"

#SIM.outputConfig.forceDD4HEP = True
#SIM.macroFile ='vis.mac'


