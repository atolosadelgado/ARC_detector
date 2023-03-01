//----------------------------------
//         ARC detector v0
//----------------------------------

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include <XML/Helper.h>

using namespace dd4hep;
using namespace dd4hep::rec;
using dd4hep::SubtractionSolid;

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

// local stuff inside anonymous namespace to avoid collisions
namespace
{
  // TODO: move parameters per cell to struct
  // store parameters by name
  std::map<std::string, double> cell_parameters_m;

  /// tokenize string by space as delimiter
  void mytokenizer(std::string &istring, std::vector<std::string> &tokens)
  {
    std::stringstream myline_ss(istring);
    std::string intermediate;
    while (getline(myline_ss, intermediate, ' '))
      tokens.push_back(intermediate);
  }

  // function to fill map with parameters cell_parameters_m
  void fill_cell_parameters_m()
  {
    // avoid calling this function twice
    if (cell_parameters_m.size())
      return;

    // hardcoded, to be developed later
    std::ifstream ifile("RadiatorCell_FinaOptimisation.txt");

    // prepare some counters for later sanity check
    int Curvature_counter(0);
    int XPosition_counter(0);
    int ZPosition_counter(0);
    int DetPosition_counter(0);
    int DetTilt_counter(0);
    int barrel_unique_cells(0);
    int endcap_unique_cells(0);

    while (ifile.good())
    {
      // read one line and tokenize by delimiter
      std::string myline("");
      std::getline(ifile, myline);
      // skip if empty line of line start with #
      if ((0 == myline.size()) || (myline.size() && '#' == myline[0]))
        continue;

      std::vector<std::string> tokens;
      mytokenizer(myline, tokens);
      // skip if not 2 elements are provided
      if (2 != tokens.size())
        continue;

      std::string &parname = tokens.at(0);
      cell_parameters_m[parname] = atof(tokens[1].c_str());

      // increase corresponding parameter counter
      // and calibrate parameter according to Martin units
      if (std::string::npos != parname.find("Curvature"))
      {
        ++Curvature_counter;
        cell_parameters_m[parname] = cell_parameters_m[parname] * m;
      }
      else if (std::string::npos != parname.find("XPosition"))
      {
        ++XPosition_counter;
        cell_parameters_m[parname] = cell_parameters_m[parname] * m;
      }
      else if (std::string::npos != parname.find("ZPosition"))
      {
        ++ZPosition_counter;
        cell_parameters_m[parname] = cell_parameters_m[parname] * m;
      }
      else if (std::string::npos != parname.find("DetPosition"))
      {
        ++DetPosition_counter;
        cell_parameters_m[parname] = cell_parameters_m[parname] * m;
      }
      else if (std::string::npos != parname.find("DetTilt"))
      {
        ++DetTilt_counter;
        cell_parameters_m[parname] = cell_parameters_m[parname] * rad;
      }
      if (std::string::npos != parname.find("EndCapRadiator"))
        ++endcap_unique_cells;
      else if (std::string::npos != parname.find("Radiator"))
        ++barrel_unique_cells;
    }
    ifile.close();

    // normalize to the number of parameters per cell
    endcap_unique_cells /= 5;
    barrel_unique_cells /= 5;

    // check if number of parameters is ok, if not, throw exception
    if (23 != endcap_unique_cells)
      throw std::runtime_error("Number of endcap cells different from expected (23)");
    if (18 != barrel_unique_cells)
      throw std::runtime_error("Number of barrel cells different from expected (18)");
    if (0 != Curvature_counter - endcap_unique_cells - barrel_unique_cells)
      throw std::runtime_error("Number of Curvature parameters different from expected (23+18)");
    if (0 != XPosition_counter - endcap_unique_cells - barrel_unique_cells)
      throw std::runtime_error("Number of XPosition parameters different from expected (23+18)");
    if (0 != ZPosition_counter - endcap_unique_cells - barrel_unique_cells)
      throw std::runtime_error("Number of ZPosition parameters different from expected (23+18)");
    if (0 != DetPosition_counter - endcap_unique_cells - barrel_unique_cells)
      throw std::runtime_error("Number of DetPosition parameters different from expected (23+18)");
    if (0 != DetTilt_counter - endcap_unique_cells - barrel_unique_cells)
      throw std::runtime_error("Number of DetTilt parameters different from expected (23+18)");

    return;
  } // end void fill_cell_parameters_m()

} // end anonymous namespace

// create barrel as sum of single cells
// next step is to place mirro+detector in a cylindral shape gas volume
static Ref_t create_barrel_cell(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  xml::Component dims = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  // read Martin file and store parameters by name in the map
  fill_cell_parameters_m();

  // mother volume corresponds to the world
  Volume motherVol = desc.pickMotherVolume(det);

  // Vessel, cylindral
  double vessel_outer_r = 210 * cm;
  double vessel_inner_r = 190 * cm;
  double vessel_length = 440 * cm;
  double vessel_wall_thickness = 1.0 * cm;
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");

  // Cooling,
  double cooling_radial_thickness = 1 * cm;

  // Cell parameters
  /// Cell is intersection of hexagonal pyramid and the cylinder
  double hexagon_side_length = 14.815 * cm;
  /// Distance in x-direction
  double zstep = 2 * hexagon_side_length;
  /// Distance in phi angle between cells
  /// since cells are regular hexagons, this distance matches
  /// the angle that one cell covers
  double phistep = 13.333 * deg;
  /// number of repetition of unique cells around the barrel
  int phinmax = 27; // 27;

  // Mirror parameters
  double thickness_sphere(10 * mm);
  double mirror_z_origin_Martin = vessel_outer_r - vessel_wall_thickness - 37 * cm;

  // Light sensor parameters
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_origin_Martin = vessel_inner_r + vessel_wall_thickness + 0.5 * cooling_radial_thickness;

  // Build cylinder for gas.
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length);

  // Use regular polyhedra for endcaps cells
  // PolyhedraRegular cellS("aa",6,0,4*cm);

  // Use pyramid for barrel cells
  std::vector<double> zplanes = {0 * cm, vessel_outer_r - vessel_wall_thickness};
  std::vector<double> rs = {0 * cm, hexagon_side_length};
  /// Hexagonal pyramid
  Polyhedra shape("mypyramid", 6, 30 * deg, 360 * deg, zplanes, rs);
  /// rotation of 90deg around Y axis, to align Z axis of pyramid with X axis of cylinder
  Transform3D pyramidTr(RotationZYX(0, -90. * deg, 0. * deg), Translation3D(0, 0, 0));

  // Build sensor shape
  Box sensorSol(sensor_sidex / 2, sensor_sidey / 2, sensor_thickness / 2);
  Volume sensorVol(detName + "_sensor", sensorSol, desc.material("Aluminum"));

  // Build the mirror for ncell=1..18
  std::vector<int> ncell_vector = {-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18,
                                   1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
  // std::vector<int> ncell_vector = { /*-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13, -14,-15,-16,-17,-18,*/
  //                                 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
  //                                 };
  // ncell_vector = {1};
  for (auto ncell : ncell_vector)
  {
    // The following line skips even number cells
    // if (!(ncell % 2))
    //   continue;

    // The following line skips odd number cells
    // if ((ncell % 2))
    // continue;

    /// cell shape. Coordinate system still the same as cylinder!
    Solid cellS = IntersectionSolid(gasvolSolid, shape, pyramidTr);
    Volume cellVol(detName + "_cell" + std::to_string(ncell), cellS, desc.material("C4F10_PFRICH"));
    cellVol.setVisAttributes(desc.visAttributes("gas_vis"));

    // there is no cell number 0, and cell number 1 do not need to be reflected
    if (0 == ncell || -1 == ncell)
      continue;

    // cells at z
    bool reflect_parameters = false;
    if (0 > ncell)
    {
      ncell *= -1;
      reflect_parameters = true;
    }

    // initialize parameters for creating the mirror
    double center_of_sphere_x(-999.);
    double center_of_sphere_z(-999.);
    double radius_of_sphere(-999.);

    double center_of_sensor_x(-999.);
    double angle_of_sensor(-999.);

    // convert Roger nomenclature (one cell number) to Martin nomenclature (row and col numbers)
    int name_col = ncell / 2;
    int name_row = ncell % 2 ? 1 : 2;
    // retrieve stored parameters
    {
      std::string name_col_s = std::to_string(name_col);
      std::string name_row_s = std::to_string(name_row);
      radius_of_sphere = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_Curvature"];
      center_of_sphere_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_XPosition"];
      double zposition = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_ZPosition"];

      center_of_sensor_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetPosition"];
      angle_of_sensor = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetTilt"];
      center_of_sphere_z = mirror_z_origin_Martin + zposition;

      // check if parameters are ok
      if (-999. == center_of_sphere_x)
        throw std::runtime_error("Ilegal parameters: center_of_sphere_x not provided");
      if (-999. == center_of_sphere_z)
        throw std::runtime_error("Ilegal parameters: center_of_sphere_z not provided");
      if (-999. == radius_of_sphere)
        throw std::runtime_error("Ilegal parameters: radius_of_sphere not provided");
      if (radius_of_sphere <= thickness_sphere)
        throw std::runtime_error("Ilegal parameters: radius_of_sphere <= thickness_sphere");

      if (-999. == center_of_sensor_x)
        throw std::runtime_error("Ilegal parameters: center_of_sensor_x not provided");
      if (-999. == angle_of_sensor)
        throw std::runtime_error("Ilegal parameters: angle_of_sensor not provided");
    }

    if (reflect_parameters)
    {
      center_of_sphere_x *= -1.0;
      center_of_sensor_x *= -1.0;
      angle_of_sensor *= -1.0;
    }

    // create the semi-sphere that will result in the mirror
    Sphere mirrorShapeFull(radius_of_sphere - thickness_sphere,
                           radius_of_sphere,
                           0.,
                           3.14 / 2);
    /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
    Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(center_of_sphere_x, 0, center_of_sphere_z));

    // TODO: cell 18 corresponds to half a pyramid, currently is full pyramid
    /// Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
    Solid mirrorSol = IntersectionSolid(shape, mirrorShapeFull, mirrorTr);
    Volume mirrorVol(detName + "_mirror" + std::to_string(ncell) + "z" + std::to_string(reflect_parameters), mirrorSol, desc.material("Aluminum"));
    mirrorVol.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell)));
    cellVol.placeVolume(mirrorVol, pyramidTr);

    // Place detector in cell
    Transform3D sensorTr(RotationZYX(0, 90 * deg - angle_of_sensor, 0), Translation3D(-sensor_z_origin_Martin, 0, center_of_sensor_x));
    cellVol.placeVolume(sensorVol, sensorTr);

    // position of mirror in cylinder coordinate system
    double mirror_abs_pos_z = name_col * zstep - 0.5 * zstep * (2 == name_row);
    if (reflect_parameters)
      mirror_abs_pos_z *= -1.0;

    // row 2 is shifted half step size
    double phi_offset = 0 + 0.5 * phistep * (2 == name_row);

    for (int phin = 0; phin < phinmax; ++phin)
    {
      PlacedVolume cellPV = motherVol.placeVolume(cellVol, RotationZ(phistep * phin + phi_offset) * Translation3D(0, 0, mirror_abs_pos_z));
      cellPV.addPhysVolID("system", detID).addPhysVolID("module", 17 * phin + name_col);
      // create mirrors as separate detectors, so properties can be adjusted lated!
      det.setPlacement(cellPV);
    }
  }

  return det;
}
DECLARE_DETELEMENT(ARCBARREL_T, create_barrel_cell)

// create endcap as sum of single cells
// next step is to place mirro+detector in a cylindral shape gas volume
static Ref_t create_endcap_cell(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  xml::Component dims = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  // read Martin file and store parameters by name in the map
  fill_cell_parameters_m();

  // mother volume corresponds to the world
  Volume motherVol = desc.pickMotherVolume(det);

  // Vessel, endcap
  double vessel_outer_r = 190 * cm;
  double vessel_inner_r = 30.2 * cm;
  double vessel_length = 20 * cm;
  double vessel_wall_thickness = 1.0 * cm;
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");

  // Cooling,
  double cooling_thickness = 1 * cm;

  // Cell parameters
  /// Cell is an hexagonal prysm
  double hexagon_side_length = 14.815 * cm;
  double hexagon_apothem = hexagon_side_length * cos(30 * deg);
  /// each cell corresponds to one object of the following class
  /// which gathers all the important parameters
  /// this info can be moved later to compact file
  struct mycell_t
  {
    /// Martin number for row
    int row = {0};
    /// Martin number for column
    int col = {0};
    /// Roger ID number
    int RID = {-1};
    /// x position of the cell
    double x = {0.};
    /// y position of the cell
    double y = {0.};
    /// if reflected
    bool isReflected = {false};
  };

  /// vector with cell geometric parameters
  std::vector<mycell_t> mycell_v(21);
  {
    double hx_u = hexagon_apothem;
    double hx_x = hexagon_side_length;
    mycell_v[0] = {1, 2, 2, 0, 4 * hx_u};
    mycell_v[1] = {1, 3, 4, 0, 6 * hx_u};
    mycell_v[2] = {1, 4, 7, 0, 8 * hx_u};
    mycell_v[3] = {1, 5, 10, 0, 10 * hx_u};
    mycell_v[4] = {1, 6, 14, 0, 12 * hx_u};
    mycell_v[5] = {1, 7, 18, 0, 14 * hx_u};
    mycell_v[6] = {2, 2, 1, -1.5 * hx_x, 3 * hx_u};
    mycell_v[7] = {2, 3, 3, -1.5 * hx_x, 5 * hx_u, true};
    mycell_v[8] = {2, 4, 6, -1.5 * hx_x, 7 * hx_u, true};
    mycell_v[9] = {2, 5, 9, -1.5 * hx_x, 9 * hx_u, true};
    mycell_v[10] = {2, 6, 13, -1.5 * hx_x, 11 * hx_u, true};
    mycell_v[11] = {2, 7, 17, -1.5 * hx_x, 13 * hx_u, true};
    mycell_v[12] = {3, 3, 5, -3.0 * hx_x, 6 * hx_u};
    mycell_v[13] = {3, 4, 8, -3.0 * hx_x, 8 * hx_u, true};
    mycell_v[14] = {3, 5, 12, -3.0 * hx_x, 10 * hx_u, true};
    mycell_v[15] = {3, 6, 16, -3.0 * hx_x, 12 * hx_u, true};
    mycell_v[16] = {3, 7, 21, -3.0 * hx_x, 14 * hx_u, true};
    mycell_v[17] = {4, 5, 11, -4.5 * hx_x, 9 * hx_u};
    mycell_v[18] = {4, 6, 15, -4.5 * hx_x, 11 * hx_u, true};
    mycell_v[19] = {4, 7, 20, -4.5 * hx_x, 13 * hx_u, true};
    mycell_v[20] = {5, 6, 19, -6.0 * hx_x, 12 * hx_u};
  }

  /// Distance in phi angle between cells
  /// since cells are regular hexagons, this distance matches
  /// the angle that one cell covers
  double phistep = 60 * deg;
  /// number of repetition of unique cells around the endcap
  int phinmax = 6; // 27;

  // Mirror parameters
  double thickness_sphere(10 * mm);
  double mirror_z_origin_Martin = vessel_length / 2. - vessel_wall_thickness - 37 * cm;

  // Light sensor parameters
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_origin_Martin = -vessel_length / 2. + vessel_wall_thickness + 0.5 * cooling_thickness;

  // Build cylinder for gas.
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length);

  // Use regular polyhedra for endcaps cells
  PolyhedraRegular cellS(6, 0, 0., hexagon_apothem, vessel_length);

  // Build sensor shape
  // Box sensorSol(sensor_sidex/2, sensor_sidey/2, sensor_thickness/2);
  // Volume sensorVol(detName + "_sensor", sensorSol, desc.material("Aluminum"));

  // Build the mirror for ncell=1..21
  // auto ncell = mycell_v[0];
  for (auto & ncell : mycell_v)
  { 
    if( -1 == ncell.RID )
      continue;
    for(int phin =0;phin<1; phin++)
    {  
    
    std::string volname = detName + "_cell" + std::to_string(ncell.RID);
    volname +="_phin" + std::to_string(phin);
    Volume cellV(volname, cellS, desc.material("C4F10_PFRICH"));


        // cellV.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell.RID)));

    // The following line skips even number cells
    // if (!(ncell % 2))
    //   continue;

    // The following line skips odd number cells
    // if ((ncell % 2))
    // continue;

    // // initialize parameters for creating the mirror
    double center_of_sphere_x(-999.);
    double center_of_sphere_z(-999.);
    double radius_of_sphere(-999.);

    // double center_of_sensor_x(-999.);
    // double angle_of_sensor(-999.);

    // convert Roger nomenclature (one cell number) to Martin nomenclature (row and col numbers)
    int name_col = ncell.col;
    int name_row = ncell.row;
    //retrieve stored parameters
    {
      std::string name_col_s = std::to_string(name_col);
      std::string name_row_s = std::to_string(name_row);
      radius_of_sphere = cell_parameters_m.at("EndCapRadiator_c" + name_col_s + "_r" + name_row_s + "_Curvature");
      center_of_sphere_x = cell_parameters_m.at("EndCapRadiator_c" + name_col_s + "_r" + name_row_s + "_XPosition");
      double zposition = cell_parameters_m.at("EndCapRadiator_c" + name_col_s + "_r" + name_row_s + "_ZPosition");
      center_of_sphere_z = mirror_z_origin_Martin + zposition;

    //   center_of_sensor_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetPosition"];
    //   angle_of_sensor = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetTilt"];

      // check if parameters are ok
      if ( -999. == center_of_sphere_x)
        throw std::runtime_error("Ilegal parameters: center_of_sphere_x not provided");
      if ( -999. == center_of_sphere_z)
        throw std::runtime_error("Ilegal parameters: center_of_sphere_z not provided");
      if ( -999. == radius_of_sphere)
        throw std::runtime_error("Ilegal parameters: radius_of_sphere not provided");
      if (radius_of_sphere <= thickness_sphere)
        throw std::runtime_error(Form("Ilegal parameters cell %d: %g <= %g", ncell.RID,radius_of_sphere, thickness_sphere ));

    //   if ( -999. == center_of_sensor_x)
    //     throw std::runtime_error("Ilegal parameters: center_of_sensor_x not provided");
    //   if ( -999. == angle_of_sensor)
    //     throw std::runtime_error("Ilegal parameters: angle_of_sensor not provided");

    }

    // if (reflect_parameters)
    // {
    //   center_of_sphere_x *= -1.0;
    //   center_of_sensor_x *= -1.0;
    //   angle_of_sensor *= -1.0;
    // }

    // create the semi-sphere that will result in the mirror
    Sphere mirrorShapeFull(radius_of_sphere - thickness_sphere,
                           radius_of_sphere,
                           0.,
                           3.14 / 2);
    /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
    double alpha = atan( ncell.y / ncell.x) * rad ;
    if( 0 > alpha)
      alpha += 180*deg;
    double dx = center_of_sphere_x*cos(alpha);
    double dy = center_of_sphere_x*sin(alpha);
    std::cout << ncell.RID << '\t' << center_of_sphere_x << '\t' << alpha/rad*deg << std::endl;
    
    Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(dx, dy, center_of_sphere_z));

    // TODO: cell 18 corresponds to half a pyramid, currently is full pyramid
    /// Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
    Solid mirrorSol = IntersectionSolid(cellS, mirrorShapeFull, mirrorTr);
    Volume mirrorVol(detName + "_mirror" + std::to_string(ncell.RID) + "z" + std::to_string(ncell.isReflected), mirrorSol, desc.material("Aluminum"));
    mirrorVol.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell.RID)));
    cellV.placeVolume(mirrorVol);

    // // Place detector in cell
    // Transform3D sensorTr(RotationZYX(0, 90*deg-angle_of_sensor, 0), Translation3D(-sensor_z_origin_Martin, 0, center_of_sensor_x));
    // cellVol.placeVolume(sensorVol, sensorTr);

    // // position of mirror in cylinder coordinate system
    // double mirror_abs_pos_z = name_col * zstep - 0.5 * zstep * (2 == name_row);
    // if (reflect_parameters)
    //   mirror_abs_pos_z *= -1.0;

    // // row 2 is shifted half step size
    // double phi_offset = 0 + 0.5 * phistep * (2 == name_row);

    // for (int phin = 0; phin < phinmax; ++phin)
    // {
    //   PlacedVolume cellPV = motherVol.placeVolume(cellVol, RotationZ(phistep * phin + phi_offset) * Translation3D(0, 0, mirror_abs_pos_z));
    //   cellPV.addPhysVolID("system", detID).addPhysVolID("module", 17 * phin + name_col);
    //   // create mirrors as separate detectors, so properties can be adjusted lated!
    //   det.setPlacement(cellPV);
    // }

    cellV.setVisAttributes(desc.visAttributes("gas_vis"));
    PlacedVolume cellPV = motherVol.placeVolume(cellV, RotationZ(phistep * phin ) * Translation3D(ncell.x, ncell.y, 0));
    cellPV.addPhysVolID("system", detID).addPhysVolID("module", ncell.RID);
    // create mirrors as separate detectors, so properties can be adjusted lated!
    det.setPlacement(cellPV);
    // if (ncell.isReflected)
    // {
    //   motherVol.placeVolume(cellV, RotationZ(phistep * phin ) * Translation3D(-ncell.x, ncell.y, 0));
    // }
  } //-- end loop for sector
  } //-- end loop for endcap
  

  return det;
}
DECLARE_DETELEMENT(ARCENDCAP_T, create_endcap_cell)

/// This fcn just build the individual cell volume, without elements
static Ref_t create_endcap_cell_volumes(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  xml::Component dims = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  // read Martin file and store parameters by name in the map
  fill_cell_parameters_m();

  // mother volume corresponds to the world
  Volume motherVol = desc.pickMotherVolume(det);

  // Vessel, endcap
  double vessel_outer_r = 190 * cm;
  double vessel_inner_r = 30.2 * cm;
  double vessel_length = 20 * cm;
  double vessel_wall_thickness = 1.0 * cm;
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");

  // Cooling,
  double cooling_thickness = 1 * cm;

  // Cell parameters
  /// Cell is an hexagonal prysm
  double hexagon_side_length = 14.815 * cm;
  double hexagon_apothem = hexagon_side_length * cos(30 * deg);
  /// each cell corresponds to one object of the following class
  /// which gathers all the important parameters
  /// this info can be moved later to compact file
  struct mycell_t
  {
    /// Martin number for row
    int row = {0};
    /// Martin number for column
    int col = {0};
    /// Roger ID number
    int RID = {0};
    /// x position of the cell
    double x = {0.};
    /// y position of the cell
    double y = {0.};
    /// if reflected
    bool isReflected = {false};
  };

  /// vector with cell geometric parameters
  std::vector<mycell_t> mycell_v(21);
  {
    double hx_u = hexagon_apothem;
    double hx_x = hexagon_side_length;
    mycell_v[0] = {1, 1, 2, 0, 4 * hx_u};
    mycell_v[1] = {1, 2, 4, 0, 6 * hx_u};
    mycell_v[2] = {1, 3, 7, 0, 8 * hx_u};
    mycell_v[3] = {1, 4, 10, 0, 10 * hx_u};
    mycell_v[4] = {1, 5, 14, 0, 12 * hx_u};
    mycell_v[5] = {1, 6, 18, 0, 14 * hx_u};
    mycell_v[6] = {2, 1, 1, -1.5 * hx_x, 3 * hx_u};
    mycell_v[7] = {2, 2, 3, -1.5 * hx_x, 5 * hx_u, true};
    mycell_v[8] = {2, 3, 6, -1.5 * hx_x, 7 * hx_u, true};
    mycell_v[9] = {2, 4, 9, -1.5 * hx_x, 9 * hx_u, true};
    mycell_v[10] = {2, 5, 13, -1.5 * hx_x, 11 * hx_u, true};
    mycell_v[11] = {2, 6, 17, -1.5 * hx_x, 13 * hx_u, true};
    mycell_v[12] = {3, 1, 5, -3.0 * hx_x, 6 * hx_u};
    mycell_v[13] = {3, 2, 8, -3.0 * hx_x, 8 * hx_u, true};
    mycell_v[14] = {3, 3, 12, -3.0 * hx_x, 10 * hx_u, true};
    mycell_v[15] = {3, 4, 16, -3.0 * hx_x, 12 * hx_u, true};
    mycell_v[16] = {3, 5, 21, -3.0 * hx_x, 14 * hx_u, true};
    mycell_v[17] = {4, 1, 11, -4.5 * hx_x, 9 * hx_u};
    mycell_v[18] = {4, 2, 15, -4.5 * hx_x, 11 * hx_u, true};
    mycell_v[19] = {4, 3, 20, -4.5 * hx_x, 13 * hx_u, true};
    mycell_v[20] = {5, 1, 19, -6.0 * hx_x, 12 * hx_u};
  }

  /// Distance in phi angle between cells
  /// since cells are regular hexagons, this distance matches
  /// the angle that one cell covers
  double phistep = 60 * deg;
  /// number of repetition of unique cells around the endcap
  int phinmax = 6; // 27;

  // Mirror parameters
  double thickness_sphere(10 * mm);
  double mirror_z_origin_Martin = vessel_length / 2. - vessel_wall_thickness - 37 * cm;

  // Light sensor parameters
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_origin_Martin = -vessel_length / 2. + vessel_wall_thickness + 0.5 * cooling_thickness;

  // Build cylinder for gas.
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length);

  // Use regular polyhedra for endcaps cells
  PolyhedraRegular cellS(6, 0, 0., hexagon_apothem, vessel_length);

  // Build the mirror for ncell=1..21

  for (auto ncell : mycell_v)
  { 
    for(int phin =0;phin<6; phin++)
    {  
    
    std::string volname = detName + "_cell" + std::to_string(ncell.RID);
    volname +="_phin" + std::to_string(phin);
    Volume cellV(volname, cellS, desc.material("Aluminum"));
    cellV.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell.RID)));
    Transform3D cellTr(RotationZ(60*phin*deg), Translation3D(ncell.x, ncell.y, 0));
    PlacedVolume cellPV = motherVol.placeVolume(cellV, RotationZ(phistep * phin ) * Translation3D(ncell.x, ncell.y, 0));
    cellPV.addPhysVolID("system", detID).addPhysVolID("module", ncell.RID);
    // create mirrors as separate detectors, so properties can be adjusted lated!
    det.setPlacement(cellPV);
    if (ncell.isReflected)
    {
      Transform3D cellTrReflected(RotationZ(60*phin*deg), Translation3D(-ncell.x, ncell.y, 0));
      motherVol.placeVolume(cellV, RotationZ(phistep * phin ) * Translation3D(-ncell.x, ncell.y, 0));
    }

  } //-- end loop for sector
  } //-- end loop for endcap
  

  return det;
}

// create minimal working example
static Ref_t createDetector(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  xml::Component dims = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  ///----------->>> constant attributes
  // - vessel
  double cell_x = dims.attr<double>(_Unicode(cell_x));
  double cell_y = dims.attr<double>(_Unicode(cell_y));
  double cell_z = dims.attr<double>(_Unicode(cell_z));
  double cell_wall_thickness = dims.attr<double>(_Unicode(cell_wall_thickness));

  auto vesselMat = desc.material(detElem.attr<std::string>(_Unicode(material)));
  auto gasvolMat = desc.material(detElem.attr<std::string>(_Unicode(gas)));
  auto vesselVis = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_vessel)));
  auto gasvolVis = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_gas)));

  // - mirror
  auto mirrorElem = detElem.child(_Unicode(mirror)).child(_Unicode(module));
  auto mirrorVis = desc.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
  auto mirrorMat = desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
  auto mirrorR = mirrorElem.attr<double>(_Unicode(radius));
  auto mirrorThickness = mirrorElem.attr<double>(_Unicode(thickness));
  auto mirrorSurf = surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));

  // - sensor module
  auto sensorElem = detElem.child(_Unicode(sensors)).child(_Unicode(module));
  auto sensorVis = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  double sensorX = sensorElem.attr<double>(_Unicode(sensorX));
  double sensorY = sensorElem.attr<double>(_Unicode(sensorY));
  double sensorThickness = sensorElem.attr<double>(_Unicode(thickness));
  auto sensorSurf = surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));
  auto sensorMat = desc.material(sensorElem.attr<std::string>(_Unicode(material)));

  // - aerogel
  auto aerogelElem = detElem.child(_Unicode(aerogel)).child(_Unicode(module));
  auto aerogelVis = desc.visAttributes(aerogelElem.attr<std::string>(_Unicode(vis)));
  auto aerogelMat = desc.material(aerogelElem.attr<std::string>(_Unicode(material)));
  double aerogel_thickness = aerogelElem.attr<double>(_Unicode(thickness));

  // - cooling
  auto coolingElem = detElem.child(_Unicode(cooling)).child(_Unicode(module));
  auto coolingVis = desc.visAttributes(coolingElem.attr<std::string>(_Unicode(vis)));
  auto coolingMat = desc.material(coolingElem.attr<std::string>(_Unicode(material)));
  double cooling_thickness = coolingElem.attr<double>(_Unicode(thickness));

  ///----------->>> Define vessel and gas volumes
  /* - `vessel`: aluminum enclosure, mother volume
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   as children of `gasvol`. Sensor is placed inside Cooling.
   * vessel (cubic) -> Gasvol (cubic) -> Cooling (cubic) -> Sensor CCD (cubic)
   *                                 \-> Mirror (sphere intersection with gasvol)
   *                                 \-> Aerogel
   */

  // Vessel
  Box vesselSolid(cell_x / 2.,
                  cell_y / 2.,
                  cell_z / 2.);
  Volume vesselVol(detName, vesselSolid, vesselMat);
  vesselVol.setVisAttributes(vesselVis);

  // Gas
  // Thickness of gas volume (z-direction) if we ignore the mirror
  double gasThickness = cell_z - 2 * cell_wall_thickness;
  Box gasvolSolid(cell_x / 2. - cell_wall_thickness,
                  cell_y / 2. - cell_wall_thickness,
                  gasThickness / 2.);
  Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
  gasvolVol.setVisAttributes(gasvolVis);
  /* PlacedVolume gasvolPV = */ vesselVol.placeVolume(gasvolVol, Position(0, 0, 0));

  ///----------->>> Aerogel (+sensor)
  {
    Box coolingSolid(cell_x / 2. - cell_wall_thickness,
                     cell_y / 2. - cell_wall_thickness,
                     cooling_thickness / 2.);
    Volume coolingVol(detName + "cooling", coolingSolid, coolingMat);
    coolingVol.setVisAttributes(coolingVis);

    // place cooling volume
    double coolingCentre = -gasThickness / 2. + cooling_thickness / 2.;
    PlacedVolume coolingPV = gasvolVol.placeVolume(coolingVol, Position(0, 0, coolingCentre));
    coolingPV.addPhysVolID("module", 63);

    ///----------->>> Sensor
    {
      Box sensorShape(sensorX / 2.,
                      sensorY / 2.,
                      sensorThickness / 2.);
      Volume sensorVol(detName + "_sensor", sensorShape, sensorMat);

      sensorVol.setVisAttributes(sensorVis);
      sensorVol.setSensitiveDetector(sens);
      double sensorCentre = cooling_thickness / 2. - sensorThickness / 2.;
      PlacedVolume sensorPV = coolingVol.placeVolume(sensorVol, Position(0, 0, sensorCentre));
      sensorPV.addPhysVolID("module", 127);

      // // Make sensor sensitive + define optical properties
      // DetElement sensorDE(aerogelDE, "ARC_sensor", 127);
      // sensorDE.setPlacement(sensorPV);
      // SkinSurface sensorSkin(desc, sensorDE, "sensor_optical_surface", sensorSurf, sensorVol); // FIXME: 3rd arg needs `imod`?
      // sensorSkin.isValid();
    }
  }
  ///----------->>> Aerogel (+sensor)
  {
    Box aerogelSolid(cell_x / 2. - cell_wall_thickness,
                     cell_y / 2. - cell_wall_thickness,
                     aerogel_thickness / 2.);
    Volume aerogelVol(detName + "_aerogel", aerogelSolid, aerogelMat);
    aerogelVol.setVisAttributes(aerogelVis);

    // place aerogel volume
    // z-position of gas volume
    double aerogelCentre = cooling_thickness - gasThickness / 2. + aerogel_thickness / 2.;
    /*PlacedVolume aerogelPV = */ gasvolVol.placeVolume(aerogelVol, Position(0, 0, aerogelCentre));
  }
  ///----------->>> Mirror
  {
    // define "mirrorVolFull" as a hollow sphere of Aluminium
    Sphere mirrorShapeFull(mirrorR - mirrorThickness,
                           mirrorR,
                           0.,
                           3.14 / 2);

    // 3D transformation of mirrorVolFull in order to place it inside the gas volume
    Transform3D mirrorTr(RotationZYX(0., 0, 0.), Translation3D(0, 0, gasThickness / 2. - mirrorR));

    // Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
    Solid mirrorSol = IntersectionSolid(gasvolSolid, mirrorShapeFull, mirrorTr);
    Volume mirrorVol(detName + "_mirror", mirrorSol, mirrorMat);
    mirrorVol.setVisAttributes(mirrorVis);
    PlacedVolume mirrorPV = gasvolVol.placeVolume(mirrorVol);
    mirrorPV.addPhysVolID("module", 3);
    DetElement mirrorDE(det, "ARC_mirror", 3);
    mirrorDE.setPlacement(mirrorPV);
    SkinSurface mirrorSkin(desc, mirrorDE, "mirror_optical_surface", mirrorSurf, mirrorVol); // FIXME: 3rd arg needs `imod`?
    mirrorSkin.isValid();
  }

  // place mother volume (vessel)
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVol, Position(0, 0, 0));
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);

  return det;
}

// clang-format off
DECLARE_DETELEMENT(ARCTYPE, createDetector)
