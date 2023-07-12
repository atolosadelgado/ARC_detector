"""
DD4hep simulation with some argument parsing
Based on M. Frank and F. Gaede runSim.py
   @author  A.Sailer
   @version 0.1

Modified with settings for RICH simulation
"""
from __future__ import absolute_import, unicode_literals
import logging
import sys
import os

from DDSim.DD4hepSimulation import DD4hepSimulation


if __name__ == "__main__":
    logging.basicConfig(
        format="%(name)-16s %(levelname)s %(message)s",
        level=logging.INFO,
        stream=sys.stdout,
        )
    logger = logging.getLogger("DDSim")

    SIM = DD4hepSimulation()

    # Default is enable visualization of tracks
    SIM.runType = "qt"
    SIM.macroFile ='vis.mac'

    # Ensure that Cerenkov and optical physics are always loaded
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

    # Allow energy depositions to 0 energy in trackers (which include optical detectors)
    SIM.filter.tracker = "edep0"

    # Some detectors are only sensitive to optical photons
    SIM.filter.filters["opticalphotons"] = dict(
        name="ParticleSelectFilter/OpticalPhotonSelector",
        parameter={"particle": "opticalphoton"},
        )
    SIM.filter.mapDetFilter["ARC_DETECTORNAME"] = "opticalphotons"

    # Use the optical tracker for the PFRICH
    SIM.action.mapActions["ARC_DETECTORNAME"] = "Geant4OpticalTrackerAction"

    # Disable user tracker particle handler, so hits can be associated to photons
    SIM.part.userParticleHandler = ""

    # Particle gun settings: pions with fixed energy and theta, varying phi
    SIM.numberOfEvents = 100
    SIM.enableGun = True
    SIM.gun.energy = "50*GeV"
    SIM.gun.particle = "pi+"
    SIM.gun.thetaMin = ".0*deg"
    SIM.gun.thetaMax = "0.1*deg"
    SIM.gun.distribution = "cos(theta)"
    SIM.gun.multiplicity = 1
    SIM.gun.position = "0 0 -20*cm"


    # Default compact file
    SIM.compactFile = "./compact/arc_v0.xml"

    # Output file (assuming CWD)
    SIM.outputFile = "arcsim.root"

    # Override with user options
    SIM.parseOptions()

    # Run the simulation
    try:
        SIM.run()
        logger.info("TEST: passed")
    except NameError as e:
        logger.fatal("TEST: failed")
        if "global name" in str(e):
            globalToSet = str(e).split("'")[1]
            logger.fatal("Unknown global variable, please add\nglobal %s\nto your steeringFile" % globalToSet)
