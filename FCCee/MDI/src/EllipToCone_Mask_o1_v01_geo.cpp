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
// Mask driver for elliptical cones, following k4geo_Mask_o1_v01: the same detector, section loop and placement
// at +z and -z, only the section parameters and the solids differ (see DDDetectors/EllipToConeHelpers.h).
#include "DDDetectors/EllipToConeHelpers.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"

#include <cmath>
#include <string>
#include <vector>

using dd4hep::Assembly;
using dd4hep::Detector;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::Ref_t;
using dd4hep::RotateX;
using dd4hep::RotateY;
using dd4hep::RotationX;
using dd4hep::RotationY;
using dd4hep::SensitiveDetector;
using dd4hep::Solid;
using dd4hep::Transform3D;
using dd4hep::Volume;

namespace units = dd4hep;

static Ref_t create_detector(Detector& description, xml_h xmlHandle, SensitiveDetector sens) {

  printout(dd4hep::DEBUG, "EllipToCone_Mask", "Creating Mask");

  // Access to the XML File
  xml_det_t xmlMask = xmlHandle;
  const std::string name = xmlMask.nameStr();

  //--------------------------------
  Assembly envelope(name + "_assembly");
  //--------------------------------

  dd4hep::DetElement tube(name, xmlMask.id());

  bool rotationX = false;

  // Parameters we have to know about
  dd4hep::xml::Component xmlParameter = xmlMask.child(_Unicode(parameter));
  const double crossingAngle = xmlParameter.attr<double>(_Unicode(crossingangle)) * 0.5; //  only half the angle

  if (xmlParameter.hasAttr(_Unicode(rotationX)))
    rotationX = xmlParameter.attr<bool>(_Unicode(rotationX));

  int counter = 0;
  for (xml_coll_t c(xmlMask, Unicode("section")); c; ++c, counter++) {

    xml_comp_t xmlSection(c);
    bool isSensitive = false;
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

    std::string ssensitive = "none";
    if (xmlSection.hasAttr(_U(sensitive))) {
      isSensitive = true;
      ssensitive = xmlSection.attr<std::string>(_U(sensitive));
      sens.setType(
          xmlSection.attr<std::string>(_U(sensitive))); // decide the type of SD (tracker / calorimeter) check for k4run
      printout(dd4hep::DEBUG, "sensitive in sens ", ssensitive);
    }

    std::stringstream pipeInfo;
    pipeInfo << std::setw(8) << zStart / units::mm << std::setw(8) << zEnd / units::mm << std::setw(8)
             << aInnerStart / units::mm << std::setw(8) << bInnerStart / units::mm << std::setw(8)
             << aInnerEnd / units::mm << std::setw(8) << bInnerEnd / units::mm << std::setw(8)
             << aOuterStart / units::mm << std::setw(8) << bOuterStart / units::mm << std::setw(8)
             << aOuterEnd / units::mm << std::setw(8) << bOuterEnd / units::mm << std::setw(8) << nSect
             << std::setw(8) << shapeType << std::setw(35) << volName << std::setw(15) << sectionMat.name()
             << std::setw(8) << ssensitive;

    printout(dd4hep::INFO, "EllipToCone_Mask", pipeInfo.str());

    // things which can be calculated immediately
    const double zHalf = fabs(zEnd - zStart) * 0.5;     // half z length of the cone
    const double zPosition = fabs(zEnd + zStart) * 0.5; // middle z position
    Material material = sectionMat;

    // this could mess up your geometry, so better check it
    if (not ECH::checkForSensibleGeometry(crossingAngle, shapeType)) {
      throw std::runtime_error(" EllipToCone_Mask_o1_v01_geo.cpp : checkForSensibleGeometry() failed ");
    }

    const double rotateAngle =
        ECH::getCurrentAngle(crossingAngle, shapeType); // for the placement at +z (better make it const now)
    const double mirrorAngle = M_PI - rotateAngle;      // for the "mirrored" placement at -z
    // the "mirroring" in fact is done by a rotation of (almost) 180 degrees around the y-axis

    switch (shapeType) {
    case ECH::kGeneralTraps: {
      // an elliptical cone on the z-axis made of nSect general trapezoids
      Transform3D transformer, transmirror;

      if (rotationX == true) {
        // absolute transformations for the placement in the world, rotate over X
        transformer = Transform3D(RotationX(rotateAngle), RotateX(Position(0, 0, zPosition), rotateAngle));
        transmirror = Transform3D(RotationX(mirrorAngle), RotateX(Position(0, 0, zPosition), mirrorAngle));
      } else {
        // absolute transformations for the placement in the world
        transformer = Transform3D(RotationY(rotateAngle), RotateY(Position(0, 0, zPosition), rotateAngle));
        transmirror = Transform3D(RotationY(mirrorAngle), RotateY(Position(0, 0, zPosition), mirrorAngle));
      }

      // solids for the tube: hollow sectors
      std::vector<Solid> tubeSolid =
          ECH::constructGeneralTrapsWall(volName + "_solid", zHalf, nSect, aInnerStart, bInnerStart, aOuterStart,
                                         bOuterStart, aInnerEnd, bInnerEnd, aOuterEnd, bOuterEnd);

      for (int i = 0; i < nSect; ++i) {
        // tube consists of the material given in the XML
        Volume tubeLog0(volName + "_" + std::to_string(i), tubeSolid[i], material);
        Volume tubeLog1(volName + "_" + std::to_string(i), tubeSolid[i], material);
        if (isSensitive) {
          tubeLog0.setSensitiveDetector(sens);
          tubeLog1.setSensitiveDetector(sens);
        }
        tubeLog0.setVisAttributes(description, xmlMask.visStr());
        tubeLog1.setVisAttributes(description, xmlMask.visStr());

        // placement of the tube in the world, both at +z and -z
        PlacedVolume placed0 = envelope.placeVolume(tubeLog0, transformer);
        PlacedVolume placed1 = envelope.placeVolume(tubeLog1, transmirror);

        placed0.addPhysVolID("side", 1);
        placed0.addPhysVolID("layer", counter);
        placed1.addPhysVolID("side", -1);
        placed1.addPhysVolID("layer", counter);
      }

    } break;

    case ECH::kTessSolids: {
      // an elliptical cone on the z-axis made of one tessellated solid
      Transform3D transformer, transmirror;

      if (rotationX == true) {
        // absolute transformations for the placement in the world, rotate over X
        transformer = Transform3D(RotationX(rotateAngle), RotateX(Position(0, 0, zPosition), rotateAngle));
        transmirror = Transform3D(RotationX(mirrorAngle), RotateX(Position(0, 0, zPosition), mirrorAngle));
      } else {
        // absolute transformations for the placement in the world
        transformer = Transform3D(RotationY(rotateAngle), RotateY(Position(0, 0, zPosition), rotateAngle));
        transmirror = Transform3D(RotationY(mirrorAngle), RotateY(Position(0, 0, zPosition), mirrorAngle));
      }

      // solid for the tube: a hollow cone
      Solid tubeSolid = ECH::constructTessSolidWall(volName + "_solid", zHalf, nSect, aInnerStart, bInnerStart,
                                                    aOuterStart, bOuterStart, aInnerEnd, bInnerEnd, aOuterEnd,
                                                    bOuterEnd);

      // tube consists of the material given in the XML
      Volume tubeLog0(volName, tubeSolid, material);
      Volume tubeLog1(volName, tubeSolid, material);
      if (isSensitive) {
        tubeLog0.setSensitiveDetector(sens);
        tubeLog1.setSensitiveDetector(sens);
      }
      tubeLog0.setVisAttributes(description, xmlMask.visStr());
      tubeLog1.setVisAttributes(description, xmlMask.visStr());

      // placement of the tube in the world, both at +z and -z
      PlacedVolume placed0 = envelope.placeVolume(tubeLog0, transformer);
      PlacedVolume placed1 = envelope.placeVolume(tubeLog1, transmirror);

      placed0.addPhysVolID("side", 1);
      placed0.addPhysVolID("layer", counter);
      placed1.addPhysVolID("side", -1);
      placed1.addPhysVolID("layer", counter);

    } break;

    default: {
      throw std::runtime_error(" EllipToCone_Mask_o1_v01_geo.cpp : fatal failure !! ??  ");
    }

    } // end switch

  } // for all xmlSections

  //--------------------------------------
  Volume mother = description.pickMotherVolume(tube);
  PlacedVolume pv(mother.placeVolume(envelope));
  pv.addPhysVolID("system", xmlMask.id()); //.addPhysVolID("side", 0 ) ;

  tube.setVisAttributes(description, xmlMask.visStr(), envelope);

  tube.setPlacement(pv);

  return tube;
}
DECLARE_DETELEMENT(EllipToCone_Mask_o1_v01, create_detector)
