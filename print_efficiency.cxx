#include "DD4hep/Detector.h"
R__LOAD_LIBRARY(libDDCore) 

void print_efficiency(){

auto lcdd = &(dd4hep::Detector::getInstance());
lcdd->fromCompact("./compact/arc_v0.xml");

const TGDMLMatrix * pp = lcdd->detector("ARC_DETECTORNAME").child("ARC_sensor").volume().material().property("APROPERTY");

std::cout << "matrix name: " << pp->GetName() << std::endl;

for(int i = 0 ; i < pp->GetRows(); ++i)
   std::cout << pp->Get(i,0)/dd4hep::eV << '\t' << pp->Get(i,1) << std::endl;


return;

}
