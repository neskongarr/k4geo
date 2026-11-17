//==========================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
//==========================================================================
// Beampipe driver for elliptical cones, following k4geo_Beampipe_o1_v01: the same envelope, section loop, nocore
// switch and hierarchy (the tube of beam material is the mother, the wall its daughter), only the section parameters
// and the solids differ (see DDDetectors/EllipToConeHelpers.h). No boolean solid is needed.
// Not taken over from k4geo_Beampipe_o1_v01: the DDRec tracking surfaces and ConicalSupportData, which are only filled
// for circular ODH::kCenter sections.

#include "DDDetectors/EllipToConeHelpers.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Printout.h"

#include "XML/Utilities.h"

#include <cmath>
#include <string>
#include <vector>

using dd4hep::Material;
using dd4hep::Position;
using dd4hep::RotateY;
using dd4hep::RotationY;
using dd4hep::Solid;
using dd4hep::Transform3D;
using dd4hep::Volume;

namespace units = dd4hep;

static dd4hep::Ref_t create_element(dd4hep::Detector& description, xml_h element, dd4hep::SensitiveDetector /*sens*/) {

  printout(dd4hep::DEBUG, "EllipToCone_Beampipe", "Creating Beampipe");

  // Access to the XML File
  xml_det_t xmlBeampipe = element;
  const std::string name = xmlBeampipe.nameStr();
  bool nocore = xmlBeampipe.nocore(false);

  dd4hep::DetElement tube(name, xmlBeampipe.id());

  // --- create an envelope volume and position it into the world ---------------------

  Volume envelope = dd4hep::xml::createPlacedEnvelope(description, element, tube);

  dd4hep::xml::setDetectorTypeFlag(element, tube);

  if (description.buildType() == dd4hep::BUILD_ENVELOPE)
    return tube;

  //-----------------------------------------------------------------------------------

  // Parameters we have to know about
  dd4hep::xml::Component xmlParameter = xmlBeampipe.child(_Unicode(parameter));
  const double crossingAngle = xmlParameter.attr<double>(_Unicode(crossingangle)) * 0.5; //  only half the angle

  for (xml_coll_t c(xmlBeampipe, Unicode("section")); c; ++c) {

    xml_comp_t xmlSection(c);

    ECH::EShapeType shapeType = ECH::getShapeType(xmlSection.attr<std::string>(_Unicode(type)));
    const double zStart = xmlSection.attr<double>(_Unicode(start));
    const double zEnd = xmlSection.attr<double>(_Unicode(end));
    // semi-axes along x (a) and y (b) of the inner and outer ellipses at the start (1) and the end (2)
    const double aInnerStart = xmlSection.attr<double>(_Unicode(aMin1));
    const double bInnerStart = xmlSection.attr<double>(_Unicode(bMin1));
    const double aInnerEnd = xmlSection.attr<double>(_Unicode(aMin2));
    const double bInnerEnd = xmlSection.attr<double>(_Unicode(bMin2));
    const double aOuterStart = xmlSection.attr<double>(_Unicode(aMax1));
    const double bOuterStart = xmlSection.attr<double>(_Unicode(bMax1));
    const double aOuterEnd = xmlSection.attr<double>(_Unicode(aMax2));
    const double bOuterEnd = xmlSection.attr<double>(_Unicode(bMax2));
    const int nSect = xmlSection.attr<int>(_Unicode(nSect)); // number of sectors (traps) or of polygon edges (tess)
    Material sectionMat = description.material(xmlSection.materialStr());
    const std::string volName = "tube_" + xmlSection.nameStr();

    std::stringstream pipeInfo;
    pipeInfo << std::setw(8) << zStart / units::mm << std::setw(8) << zEnd / units::mm << std::setw(8)
             << aInnerStart / units::mm << std::setw(8) << bInnerStart / units::mm << std::setw(8)
             << aInnerEnd / units::mm << std::setw(8) << bInnerEnd / units::mm << std::setw(8)
             << aOuterStart / units::mm << std::setw(8) << bOuterStart / units::mm << std::setw(8)
             << aOuterEnd / units::mm << std::setw(8) << bOuterEnd / units::mm << std::setw(8) << nSect
             << std::setw(8) << shapeType << std::setw(35) << volName << std::setw(15) << sectionMat.name();

    printout(dd4hep::INFO, "EllipToCone_Beampipe", pipeInfo.str());

    // things which can be calculated immediately
    const double zHalf = fabs(zEnd - zStart) * 0.5;       // half z length of the cone
    const double zPosition = fabs(zEnd + zStart) * 0.5;   // middle z position
    Material coreMaterial = description.material("beam"); // always the same
    Material wallMaterial = sectionMat;

    // this could mess up your geometry, so better check it
    if (not ECH::checkForSensibleGeometry(crossingAngle, shapeType)) {

      throw std::runtime_error(" EllipToCone_Beampipe_o1_v01_geo.cpp : checkForSensibleGeometry() failed ");
      //      return false;
    }
    const double rotateAngle =
        ECH::getCurrentAngle(crossingAngle, shapeType); // for the placement at +z (better make it const now)
    const double mirrorAngle = M_PI - rotateAngle;      // for the "mirrored" placement at -z
    // the "mirroring" in fact is done by a rotation of (almost) 180 degrees around the y-axis

    // if inner and outer semi-axes are equal, then omit the tube wall
    const bool hasWall = (aInnerStart != aOuterStart || bInnerStart != bOuterStart || aInnerEnd != aOuterEnd ||
                          bInnerEnd != bOuterEnd);

    switch (shapeType) {
    case ECH::kGeneralTraps: {
      // an elliptical cone on the z-axis made of nSect general trapezoids

      // absolute transformations for the placement in the world
      Transform3D transformer(RotationY(rotateAngle), RotateY(Position(0, 0, zPosition), rotateAngle));
      Transform3D transmirror(RotationY(mirrorAngle), RotateY(Position(0, 0, zPosition), mirrorAngle));

      // solids for the tube (including vacuum and wall): solid sectors of the outer cone
      std::vector<Solid> tubeSolid =
          ECH::constructGeneralTrapsTube(volName + "_solid", zHalf, nSect, aOuterStart, bOuterStart, aOuterEnd,
                                         bOuterEnd);

      if (hasWall) {

        // the wall solids: hollow sectors, each one inside the tube sector of the same index
        std::vector<Solid> wallSolid =
            ECH::constructGeneralTrapsWall(volName + "_wall_solid", zHalf, nSect, aInnerStart, bInnerStart,
                                           aOuterStart, bOuterStart, aInnerEnd, bInnerEnd, aOuterEnd, bOuterEnd);

        for (int i = 0; i < nSect; ++i) {
          const std::string sector = "_" + std::to_string(i);

          // tube consists of vacuum
          // place tube twice explicitely, as in k4geo_Beampipe_o1_v01
          Volume tubeLog(volName + "_0" + sector, tubeSolid[i], coreMaterial);
          Volume tubeLog2(volName + "_1" + sector, tubeSolid[i], coreMaterial);

          // the wall consists of the material given in the XML
          Volume wallLog(volName + "_wall_0" + sector, wallSolid[i], wallMaterial);
          Volume wallLog2(volName + "_wall_1" + sector, wallSolid[i], wallMaterial);

          wallLog.setVisAttributes(description, "TubeVis");
          wallLog2.setVisAttributes(description, "TubeVis");
          tubeLog.setVisAttributes(description, "VacVis");
          tubeLog2.setVisAttributes(description, "VacVis");

          if (nocore) {

            // placement of the wall only
            envelope.placeVolume(wallLog, transformer);
            envelope.placeVolume(wallLog2, transmirror);

          } else {

            // placement of the tube in the world, both at +z and -z
            envelope.placeVolume(tubeLog, transformer);
            envelope.placeVolume(tubeLog2, transmirror);

            // placement as a daughter volume of the tube, will appear in both placements of the tube
            tubeLog.placeVolume(wallLog, Transform3D());
            tubeLog2.placeVolume(wallLog2, Transform3D());
          }
        }
      }
    } break;

    case ECH::kTessSolids: {
      // an elliptical cone on the z-axis made of one tessellated solid

      // absolute transformations for the placement in the world
      Transform3D transformer(RotationY(rotateAngle), RotateY(Position(0, 0, zPosition), rotateAngle));
      Transform3D transmirror(RotationY(mirrorAngle), RotateY(Position(0, 0, zPosition), mirrorAngle));

      // solid for the tube (including vacuum and wall): a solid cone
      Solid tubeSolid =
          ECH::constructTessSolidTube(volName + "_solid", zHalf, nSect, aOuterStart, bOuterStart, aOuterEnd, bOuterEnd);

      // tube consists of vacuum
      // place tube twice explicitely, as in k4geo_Beampipe_o1_v01
      Volume tubeLog(volName + "_0", tubeSolid, coreMaterial);
      Volume tubeLog2(volName + "_1", tubeSolid, coreMaterial);

      if (hasWall) {

        // the wall solid: a hollow cone
        Solid wallSolid = ECH::constructTessSolidWall(volName + "_wall_solid", zHalf, nSect, aInnerStart, bInnerStart,
                                                      aOuterStart, bOuterStart, aInnerEnd, bInnerEnd, aOuterEnd,
                                                      bOuterEnd);

        // the wall consists of the material given in the XML
        Volume wallLog(volName + "_wall_0", wallSolid, wallMaterial);
        Volume wallLog2(volName + "_wall_1", wallSolid, wallMaterial);

        wallLog.setVisAttributes(description, "TubeVis");
        wallLog2.setVisAttributes(description, "TubeVis");
        tubeLog.setVisAttributes(description, "VacVis");
        tubeLog2.setVisAttributes(description, "VacVis");

        if (nocore) {

          // placement of the wall only
          envelope.placeVolume(wallLog, transformer);
          envelope.placeVolume(wallLog2, transmirror);

        } else {

          // placement of the tube in the world, both at +z and -z
          envelope.placeVolume(tubeLog, transformer);
          envelope.placeVolume(tubeLog2, transmirror);

          // placement as a daughter volume of the tube, will appear in both placements of the tube
          tubeLog.placeVolume(wallLog, Transform3D());
          tubeLog2.placeVolume(wallLog2, Transform3D());
        }
      }
    } break;

    default: {
      throw std::runtime_error(" EllipToCone_Beampipe_o1_v01_geo.cpp : fatal failure !! ??  ");
    }
    } // end switch
  } // for all xmlSections

  //--------------------------------------

  tube.setVisAttributes(description, xmlBeampipe.visStr(), envelope);

  return tube;
}
DECLARE_DETELEMENT(EllipToCone_Beampipe_o1_v01, create_element)
