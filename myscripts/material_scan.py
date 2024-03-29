# Code taken from Brieuc for performing material scan using Gaudi
# https://raw.githubusercontent.com/BrieucF/LAr_scripts/main/geometry/material_scan.py

from Gaudi.Configuration import *


from Configurables import ApplicationMgr
ApplicationMgr().EvtSel = 'None'
ApplicationMgr().EvtMax = 1
ApplicationMgr().OutputLevel = INFO

# DD4hep geometry service
from Configurables import GeoSvc

## parse the given xml file
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [ 'lcgeoTests/compact/ARC_standalone_o1_v01.xml' ]
geoservice.OutputLevel = INFO
ApplicationMgr().ExtSvc += [geoservice]

from Configurables import MaterialScan_2D_genericAngle
# Material scan is done from the interaction point to the end of world volume.
# In order to use other end boundary, please provide the name of a thin, e.g. cylindrical volume.
# For instance adding envelopeName="BoundaryPostCalorimetry" will perform the scan only till the end of calorimetry.
# BoundaryPostCalorimetry is defined in Detector/DetFCChhECalInclined/compact/envelopePreCalo.xml
materialservice = MaterialScan_2D_genericAngle("GeoDump")
materialservice.filename = "out_material_scan.root"
materialservice.angleMax = 0.999
materialservice.angleMin = 0
materialservice.angleBinning = 0.001
materialservice.angleDef = 'cosTheta'
materialservice.nPhi = 1
ApplicationMgr().ExtSvc += [materialservice]



