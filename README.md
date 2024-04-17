ARC RICH detector V0
=======================

Simple cell description of the ARC RICH detector.

The directory structure is the following

* README.md: this file
* CMakeLists.txt: minimal cmake macro to configure the detector construction
* compact: has material and elements with their description (including optical properties). The compact file of the detector is also here.
* src: detector consctructor, writen in c++
* myscripts: useful bash scripts and GEANT4 macro files


# Basic commands

Build with `cmake`, 
```bash
source /cvmfs/sft.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
cmake -B build -S . -D CMAKE_INSTALL_PREFIX=install
cmake --build build -- install
export LD_LIBRARY_PATH=$PWD/install/lib:$LD_LIBRARY_PATH
```

To display the geometry

```bash
geoDisplay compact/arc_v0.xml 
```

To convert the geometry into ROOT

```bash
./myscripts/dd4hep2root -c  compact/arc_v0.xml -o arc_v0.root
```

Using Geant4+Qt, particularly useful when using boolean operations and ROOT display do not work properly.

```bash
ddsim --compactFile compact/arc_barrel_v0.xml --runType qt --macroFile vis.mac --part.userParticleHandler=''
```

# Things to do after every change

Compile and install after every modification of the c++ detector constructor code.

```bash
cmake --build build -- install
```

To show materials in along a given line (in this case from origin to (0,0,-100)cm), check if the modified volume has the proper material

```bash
materialScan compact/arc_v0.xml 0 0 0 0 0 -100
```

Check overlaps, using Geant4 check
```shell
ddsim --compactFile ./compact/arc_v0.xml --runType run --part.userParticleHandler='' --macroFile myscripts/overlap.mac >> overlapDump.txt
```


# Run simulation

The script `arcsim.py` is a ddsim steering file, which enable Qt visualization, setup the cerenkov physic lists and the default input and output files. The configuration inside the steering file can be overriden if options are added after the name of the steering file (see -N 1 example below). 

Use this command to open the Geant4 Qt application,

```shell
python arcsim.py -N 1
```

then press play button on top to launch one event. The geometry and tracks of one charged pion (pink) crossing the detector, and the cerenkov photons (grey) are displayed below.


![Visualization of geometry and tracks of one charged pion (pink) crossing the detector, and the cerenkov photons (grey).](https://mattermost.web.cern.ch/files/4eyjw6q8b38tfkmh5nrf7qytsr/public?h=rLXJPn3Zsd6-0g3Q9ZaqHI9pAFglRqD_kix70QP6nXs)



More events can be simulated, but in this case only 1 will be saved into the output file `arcsim.root`. To visualize the hit pattern in the detector, we can open the root file and write the following in the command line

```cpp
EVENT->Draw("ARC_HITS.position.Y():ARC_HITS.position.X()");
```

In case of EDM4hep output format, execute the following line instead

```cpp
events->Draw("ARC_HITS.position.y:ARC_HITS.position.x");
```

A new canvas will open showing something similar to the following graph:

![Hit pattern of photons (and pion) in the detector](https://mattermost.web.cern.ch/files/11f17b5nctfkjqjsw1cmoh1rko/public?h=OQCOs1RkwC560pOU0reOnG9pJabN4rDqTu2wgqHeHNg)

Hit pattern of photons (and pion) in the detector

## Running the simulation of the full ARC detector

There is an equivalent python steering file for the whole detector, which includes the barrel and endcaps. The following command will simulate 10k protons with default energy (50 GeV), and the output format will be EDM4hep (because of the extension of the output file name).

```
python3 arcfullsim.py --runType batch -N 10000 --gun.particle "proton" --outputFile "arcsim_proton_50GeV_edm4hep.root"
```

## Interpreting the cellID

The readout id is a 64 bit field that is specified in the xml file,
```
      <id>system:5,barrel:3,cellnumber:24,x:32:-16,y:-16</id>
```

Each name define a field, which has the corresponding number of bits. System field name identifies if it is ARC barrrel or endcap (or any main sub-detector). Barrel field name identifies if it is barrel (0) or endcap (1,2) of the ARC sub-detector. The cellnumber is just a counter for each cell. Note that for the endcap the numbering repeats, but the barrel field is different. The XY field names corresponds to the sensor segmentation, which is the only sensitive component of this subdetector.


To plot the detector ID bit field
```
 EVENT->Scan("((ARC_HITS.cellID&0x1F))>>h(100,0,100)")
```

Similarly, the barrel bitfield can be shown by using the expression `((ARC_HITS.cellID>>5)&0x7)`. It can be read as jumping the first 5 bits (corresponding to the detector ID) and crop the following 3 bits, corresponding to the barrel bitfield. This value can be 0 (barrel) or (1,2) for the endcap.

To plot the cell number (it is just a consecutive number that is assigned to each cell). The expression to be plotted can read as jumping the first 8 bits (corresponding to detector ID + barrel bitfield), and crop the following 24 bits, corresponding to the cell number. The second expression is a condition, to select the events happeing in one endcap, the one at Z>0.
```
EVENT->Draw("((ARC_HITS.cellID>>8)&0xFFFFFF)>>hCellN(300,0,300)","((ARC_HITS.cellID>>5)&0x7)==1")
```

One can inspect the hit pattern in the sensor from the cellID parameters,

```
EVENT->Scan("((ARC_HITS.cellID>>32)&0xFFFF):(ARC_HITS.cellID>>48)" )
```

or draw the sensor hit pattern just the fifth event


```
 EVENT->Draw("((ARC_HITS.cellID/pow(2,32))&0xFFFF):(ARC_HITS.cellID/pow(2,48))","","colz",1,5)
```

Note 1: this primitive way of chopping the bitfield can be done using the corresponding DD4hep object.
Note 2: if EDM4hep data format is used, please replace the name of the tree `EVENT` by `events`. The `cellID` tricks will work, but other variable names may change a bit.

# Tests

A simple test can be running a simulation with 1k of pi+ at 50 GeV, as the following

```
python3 arcfullsim.py --runType batch
```

Endcaps or barrel can be tested separately as well,

```
python3 endcapsim.py --runType batch
python3 barrelsim.py --runType batch
```

Visualization of the detector can be also a test by itself,

```
geoDisplay compact/arc_barrel_v0.xml
geoDisplay compact/arc_endcap_v0.xml
geoDisplay compact/arc_full_v0.xml
```

# Migration to central detector repository k4geo

This code was transferred to the FCCee central detector repository [k4geo](https://github.com/key4hep/k4geo/tree/master/detector/PID). In order to run a DD4hep stand-alone simulation of ARC detector with the code from k4geo, please follow the following steps:

1. Download, compile and make the k4geo repository visible to the environment
```
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
git clone -b master  git@github.com:atolosadelgado/k4geo.git
cd k4geo/
cmake -B build -S . -D CMAKE_INSTALL_PREFIX=./install
cmake --build build -j 20 -- install
export LD_LIBRARY_PATH=$PWD/install/lib:$LD_LIBRARY_PATH
```

2. Move to k4geo/example directory, and you can launch a simulation of ARC with the following command

```
cd example
python3 arcfullsim.py --runType batch --compactFile ../FCCee/CLD/compact/CLD_o3_v01/CLD_o3_v01.xml --gun.particle "proton" --gun.energy "10*GeV" -N 1234 --outputFile arcsim_edm4hep.root
```

Particle type must be provided as a string, the full list of possible particles can be found here: http://fismed.ciemat.es/GAMOS/GAMOS_doc/GAMOS.5.1.0/x11519.html

3. You can further customize the particle gun with the following flags:
```
--gun.particle GUN.PARTICLE
    --gun.multiplicity GUN.MULTIPLICITY
    --gun.phiMin GUN.PHIMIN
    --gun.phiMax GUN.PHIMAX
    --gun.thetaMin GUN.THETAMIN
    --gun.thetaMax GUN.THETAMAX
    --gun.etaMin GUN.ETAMIN
    --gun.etaMax GUN.ETAMAX
    --gun.momentumMin GUN.MOMENTUMMIN
    --gun.momentumMax GUN.MOMENTUMMAX
    --gun.energy GUN.ENERGY
                          The kinetic energy for the particle gun.
    --gun.direction GUN.DIRECTION
                           direction of the particle gun, 3 vector
    --gun.distribution {uniform,cos(theta),eta,pseudorapidity,ffbar}
    --gun.isotrop GUN.ISOTROP
                           isotropic distribution for the particle gun
    --gun.position GUN.POSITION
                           position of the particle gun, 3 vector
```

If even these options are not enough, Geant4 provides the so-called GPS, which can produce primary particles in a more sophisticated way (for example, histogram sampling)

# Useful links

Documentation of DD4hep,
* PDF: https://dd4hep.web.cern.ch/dd4hep/usermanuals/DD4hepManual/DD4hepManual.pdf
* Code: https://github.com/AIDASoft/DD4hep
* Class reference: https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Detector.html
* Beginners guide: https://dd4hep.web.cern.ch/dd4hep/page/beginners-guide/
* FCC tutorial: https://hep-fcc.github.io/fcc-tutorials/full-detector-simulations/Geometry/Geometry.html#geometry-driver-modifications
* NPdet: https://eic.phy.anl.gov/npdet/index.html

Example of RICH detector implemented in DD4hep
* Proximity Focusing RICH: https://github.com/AIDASoft/DD4hep/tree/master/examples/OpticalTracker

Documentation and presentations about ARC,
* 1-July-2021 FCC week https://indico.cern.ch/event/995850/contributions/4406336/attachments/2274813/3864163/ARC-presentation.pdf
* 23-June-2022 Workshop on detector for FCCee https://indico.cern.ch/event/1165167/contributions/4927647/attachments/2468337/4233740/ARC-update.pdf
* 6-October-2022 ECFA Workshop https://indico.desy.de/event/33640/contributions/128392/attachments/77681/100501/ARC_Presentation_ECFA_Workshop_DESY_6th_October_2022.pdf
* 15-Dec-2022 ARC meeting indico: https://indico.cern.ch/event/1231098/
* 25-January-2023 FCC workshop https://indico.cern.ch/event/1176398/contributions/5208403/attachments/2581656/4452921/ARC_Presentation_FCC_Workshop_Krakow_25th_January_2023.pdf
* 28-March-2023 Discussion about ARC implementation in DD4hep: https://indico.cern.ch/event/1266428/
* 8-June-2023 FCC week https://indico.cern.ch/event/1202105/contributions/5396830/attachments/2662099/4612055/fccweek23_atd_v7.pdf
* 12-July-2023 Fullsim meeting (M. Basso) https://indico.cern.ch/event/1304946/contributions/5488302/attachments/2682619/4654996/230712_Basso_FCC-Det-FS-Meeting_ARC-reco.pdf
* 28-August-2023 Fullsim meeting https://indico.cern.ch/event/1318491/contributions/5548508/attachments/2703694/4712665/ARC_IntegrationCLD_230823.pdf
* 24-October-2023 Fullsim meeting (M. Tat) https://indico.cern.ch/event/1339694/contributions/5639529/attachments/2739715/4765514/ARC_Presentation_FCC_ImplementationDiscussion_24th_October_2023.pdf
* 24-October-2023 Fullsim meeting https://indico.cern.ch/event/1339694/contributions/5639532/attachments/2739682/4765470/ARC_updates_231024.pdf
* 1-February-2024 FCC workshop. https://indico.cern.ch/event/1307378/contributions/5729665/attachments/2789636/4867784/fccphysweek24_ARC_reco_atd_240130.pdf

FCCee fullsim documentation
* https://fcc-ee-detector-full-sim.docs.cern.ch/
