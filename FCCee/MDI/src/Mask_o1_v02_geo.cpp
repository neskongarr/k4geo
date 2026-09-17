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
#include "DDDetectors/OtherDetectorHelpers.h"
#include "DDDetectors/DetectorNesting.h"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"

#include <cmath>
#include <string>

using dd4hep::Transform3D;
using dd4hep::Position;
using dd4hep::RotationX;
using dd4hep::RotationY;
using dd4hep::RotateY;
using dd4hep::RotateX;
using dd4hep::ConeSegment;
using dd4hep::SubtractionSolid;
using dd4hep::Material;
using dd4hep::Volume;
using dd4hep::Solid;
using dd4hep::Tube;
using dd4hep::PlacedVolume;
using dd4hep::Assembly;
using dd4hep::Detector;
using dd4hep::SensitiveDetector;
using dd4hep::Ref_t;
using dd4hep::Rotation3D;
using dd4hep::RotationZ;
using dd4hep::Trd1;
using dd4hep::Trd2;

namespace units = dd4hep;


static Ref_t create_detector(Detector& description,
			     xml_h xmlHandle,
			     SensitiveDetector sens) {

  printout(dd4hep::DEBUG,"DD4hep_Mask", "Creating Mask" ) ;

  //Access to the XML File
  xml_det_t xmlMask = xmlHandle;
  const std::string name = xmlMask.nameStr();

  //--------------------------------
  Assembly envelope( name + "_assembly"  ) ;
  //--------------------------------

  dd4hep::DetElement tube(  name, xmlMask.id()  ) ;

  bool rotationX= false;

  //Parameters we have to know about
  dd4hep::xml::Component xmlParameter = xmlMask.child(_Unicode(parameter));
  const double crossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle))*0.5; //  only half the angle

  if (xmlParameter.hasAttr(_Unicode(rotationX)))
    rotationX = xmlParameter.attr< bool >(_Unicode(rotationX));

  // Optional <nest_in detector="..."/>: place this detector's sections into the
  // volume of the named host detector that geometrically contains them, instead of
  // into our own assembly envelope. Only a *detector* is named here; the mother
  // volume itself is resolved geometrically, so no internal volume name is needed.
  Volume      hostVol;      // invalid => classic behaviour, nothing changes
  Transform3D hostToWorld;  // identity
  std::string hostName;
  if ( xml_h nestIn = xmlMask.child(_Unicode(nest_in), false) )
    hostName = dd4hep::xml::Component(nestIn).attr< std::string >(_Unicode(detector));

  if ( !hostName.empty() ) {
    // Deliberately a plain map lookup: Detector::detector() falls through to a
    // DetElement path search whose failure message ("Unknown child with name") is
    // unhelpful for what is really an XML ordering mistake.
    const auto& dets = description.detectors();
    auto it = dets.find(hostName);
    if ( it == dets.end() ) {
      std::string known;
      for (const auto& d : dets) known += " " + d.first;
      dd4hep::except("MDI_Mask",
                     "+++ %s: <nest_in detector=\"%s\"/> : no such detector has been built "
                     "yet. <detector> elements are processed in XML document order, so the "
                     "host MUST appear before %s in the compact file. Built so far:%s",
                     name.c_str(), hostName.c_str(), name.c_str(), known.c_str());
    }
    dd4hep::DetElement host((*it).second);
    if ( !host.placement().isValid() )
      dd4hep::except("MDI_Mask",
                     "+++ %s: <nest_in detector=\"%s\"/> : the host exists but has no "
                     "placement, i.e. it was not completely built. Move its <detector> "
                     "element before %s.", name.c_str(), hostName.c_str(), name.c_str());
    hostVol     = host.placement().volume();
    hostToWorld = ODH::placementToWorld(host);
  }

  int counter = 0; 
  for(xml_coll_t c( xmlMask ,Unicode("section")); c; ++c, counter++) {

    xml_comp_t xmlSection( c );
    bool isSensitive = false;
    ODH::ECrossType crossType = ODH::getCrossType(xmlSection.attr< std::string >(_Unicode(type)));
    const double zStart       = xmlSection.attr< double > (_Unicode(start));
    const double zEnd         = xmlSection.attr< double > (_Unicode(end));
    const double rInnerStart  = xmlSection.attr< double > (_Unicode(rMin1));
    const double rInnerEnd    = xmlSection.attr< double > (_Unicode(rMin2));
    const double rOuterStart  = xmlSection.attr< double > (_Unicode(rMax1));
    const double rOuterEnd    = xmlSection.attr< double > (_Unicode(rMax2));
    const double thickness    = rOuterStart - rInnerStart;
    Material sectionMat  = description.material(xmlSection.materialStr());
    const std::string volName      = "tube_" + xmlSection.nameStr();

    double phi1 = 0 ;
    double phi2 = 360.0*units::degree;
    if (xmlSection.hasAttr(_U(phi1)))
      phi1 = xmlSection.attr< double > (_U(phi1));
    if (xmlSection.hasAttr(_U(phi2)))
      phi2 = xmlSection.attr< double > (_U(phi2));

    std::string ssensitive = "none";
    if (xmlSection.hasAttr(_U(sensitive))){
      isSensitive = true;
      ssensitive = xmlSection.attr< std::string > (_U(sensitive));
      sens.setType( xmlSection.attr< std::string > (_U(sensitive))  );  //decide the type of SD (tracker / calorimeter) check for k4run
      printout(dd4hep::DEBUG, "sensitive in sens ", ssensitive);
    }


    std::stringstream pipeInfo;
    pipeInfo << std::setw(8) << zStart      /units::mm
	     << std::setw(8) << zEnd        /units::mm
	     << std::setw(8) << rInnerStart /units::mm
	     << std::setw(8) << rInnerEnd   /units::mm
	     << std::setw(8) << rOuterStart /units::mm
	     << std::setw(8) << rOuterEnd   /units::mm
	     << std::setw(8) << thickness   /units::mm
	     << std::setw(8) << crossType
	     << std::setw(35) << volName
	     << std::setw(15) << sectionMat.name()
	     << std::setw(8) << phi1
	     << std::setw(8) << phi2
	     << std::setw(8) << ssensitive;

    printout(dd4hep::INFO, "DD4hep_Mask", pipeInfo.str() );

    // things which can be calculated immediately
    const double zHalf       = fabs(zEnd - zStart) * 0.5; // half z length of the cone
    const double zPosition   = fabs(zEnd + zStart) * 0.5; // middle z position
    Material material    = sectionMat;

    // this could mess up your geometry, so better check it
    if (not ODH::checkForSensibleGeometry(crossingAngle, crossType)){
      throw std::runtime_error( " Mask_o1_v01_geo.cpp : checkForSensibleGeometry() failed " ) ;
    }

    const double rotateAngle = getCurrentAngle(crossingAngle, crossType); // for the placement at +z (better make it const now)
    const double mirrorAngle = M_PI - rotateAngle; // for the "mirrored" placement at -z
    // the "mirroring" in fact is done by a rotation of (almost) 180 degrees around the y-axis

    
    switch (crossType) {
    case ODH::kCenter:
    case ODH::kUpstream:
    case ODH::kDnstream: {
      // a volume on the z-axis, on the upstream branch, or on the downstream branch
      Transform3D transformer, transmirror;

      if( rotationX == true) {
        // absolute transformations for the placement in the world, rotate over X
        transformer = Transform3D(RotationX(rotateAngle), RotateX( Position(0, 0, zPosition), rotateAngle) );
        transmirror = Transform3D(RotationX(mirrorAngle), RotateX( Position(0, 0, zPosition), mirrorAngle) );
      } else{
	// absolute transformations for the placement in the world
	transformer = Transform3D(RotationY(rotateAngle), RotateY( Position(0, 0, zPosition), rotateAngle) );
	transmirror = Transform3D(RotationY(mirrorAngle), RotateY( Position(0, 0, zPosition), mirrorAngle) );
      }
      
      // solid for the tube (including vacuum and wall): a solid cone
      ConeSegment tubeSolid( zHalf, rInnerStart, rOuterStart, rInnerEnd, rOuterEnd , phi1, phi2);

      // tube consists of vacuum
      Volume tubeLog0( volName, tubeSolid, material ) ;
      Volume tubeLog1( volName, tubeSolid, material ) ;
      if (isSensitive) {
        tubeLog0.setSensitiveDetector(sens);
        tubeLog1.setSensitiveDetector(sens);
      }
      tubeLog0.setVisAttributes(description, xmlMask.visStr() );
      tubeLog1.setVisAttributes(description, xmlMask.visStr() );

      // placement of the tube in the world, both at +z and -z
      PlacedVolume placed0 = envelope.placeVolume( tubeLog0,  transformer );
      PlacedVolume placed1 = envelope.placeVolume( tubeLog1,  transmirror );

      placed0.addPhysVolID("side",1);
      placed0.addPhysVolID("layer",counter);
      placed1.addPhysVolID("side",-1);
      placed1.addPhysVolID("layer",counter);
      
    }
      break;

    case ODH::kPunchedCenter: {
      // a cone with one or two inner holes (two tubes are punched out)

      const double rUpstreamPunch = rInnerStart; // just alias names denoting what is meant here
      const double rDnstreamPunch = rInnerEnd; // (the database entries are "abused" in this case)

      // relative transformations for the composition of the SubtractionVolumes
      Transform3D upstreamTransformer(RotationY(-crossingAngle), Position(zPosition * tan(-crossingAngle), 0, 0));
      Transform3D dnstreamTransformer(RotationY(+crossingAngle), Position(zPosition * tan(+crossingAngle), 0, 0));

      // absolute transformations for the final placement in the world (angles always equal zero and 180 deg)
      Transform3D placementTransformer(RotationY(rotateAngle), RotateY( Position(0, 0, zPosition) , rotateAngle) );
      Transform3D placementTransmirror(RotationY(mirrorAngle), RotateY( Position(0, 0, zPosition) , mirrorAngle) );

      // the main solid and the two pieces (only tubes, for the moment) which will be punched out
      ConeSegment wholeSolid(  zHalf, 0, rOuterStart, 0, rOuterEnd, phi1, phi2 );
      Solid tmpSolid0, tmpSolid1, finalSolid0, finalSolid1;

      // the punched subtraction solids can be asymmetric and therefore have to be created twice:
      // one time in the "right" way, another time in the "reverse" way, because the "mirroring"
      // rotation around the y-axis will not only exchange +z and -z, but also +x and -x

      if ( rUpstreamPunch > 1e-6 ) { // do we need a hole on the upstream branch?
	Tube upstreamPunch( 0, rUpstreamPunch, 5 * zHalf, phi1, phi2); // a bit longer
	tmpSolid0 = SubtractionSolid( wholeSolid, upstreamPunch, upstreamTransformer);
	tmpSolid1 = SubtractionSolid( wholeSolid, upstreamPunch, dnstreamTransformer); // [sic]
      } else { // dont't do anything, just pass on the unmodified shape
	tmpSolid0 = wholeSolid;
	tmpSolid1 = wholeSolid;
      }

      if (rDnstreamPunch > 1e-6 ) { // do we need a hole on the downstream branch?
	Tube dnstreamPunch( 0, rDnstreamPunch, 5 * zHalf, phi1, phi2); // a bit longer
	finalSolid0 = SubtractionSolid( tmpSolid0, dnstreamPunch, dnstreamTransformer);
	finalSolid1 = SubtractionSolid( tmpSolid1, dnstreamPunch, upstreamTransformer); // [sic]
      } else { // dont't do anything, just pass on the unmodified shape
	finalSolid0 = tmpSolid0;
	finalSolid1 = tmpSolid1;
      }

      // tube consists of vacuum (will later have two different daughters)
      Volume tubeLog0( volName + "_0", finalSolid0, material );
      Volume tubeLog1( volName + "_1", finalSolid1, material );
      if (isSensitive) {
        tubeLog0.setSensitiveDetector(sens);
        tubeLog1.setSensitiveDetector(sens);
      }
      tubeLog0.setVisAttributes(description, xmlMask.visStr() );
      tubeLog1.setVisAttributes(description, xmlMask.visStr() );

      // placement of the tube in the world, both at +z and -z
      PlacedVolume placed0 = envelope.placeVolume( tubeLog0, placementTransformer );
      PlacedVolume placed1 = envelope.placeVolume( tubeLog1, placementTransmirror );
      
      placed0.addPhysVolID("side",  1);
      placed1.addPhysVolID("side", -1);
      placed0.addPhysVolID("layer", counter);
      placed1.addPhysVolID("layer", counter);

      break;
    }

    case ODH::kPunchedUpstream:
    case ODH::kPunchedDnstream: {
      // a volume on the upstream or downstream branch with two inner holes
      // (implemented as a cone from which another tube is punched out)

      const double rCenterPunch = (crossType == ODH::kPunchedUpstream) ? (rInnerStart) : (rInnerEnd); // just alias names denoting what is meant here
      const double rOffsetPunch = (crossType == ODH::kPunchedDnstream) ? (rInnerStart) : (rInnerEnd); // (the database entries are "abused" in this case)

      // relative transformations for the composition of the SubtractionVolumes
      Transform3D punchTransformer(RotationY(-2 * rotateAngle), Position(zPosition * tan(-2 * rotateAngle), 0, 0));
      Transform3D punchTransmirror(RotationY(+2 * rotateAngle), Position(zPosition * tan(+2 * rotateAngle), 0, 0));

      // absolute transformations for the final placement in the world
      Transform3D placementTransformer(RotationY(rotateAngle), RotateY( Position(0, 0, zPosition) , rotateAngle) );
      Transform3D placementTransmirror(RotationY(mirrorAngle), RotateY( Position(0, 0, zPosition) , mirrorAngle) );

      // the main solid and the piece (only a tube, for the moment) which will be punched out
      ConeSegment wholeSolid( zHalf, rCenterPunch , rOuterStart, rCenterPunch, rOuterEnd, phi1, phi2);
      Tube punchSolid( 0, rOffsetPunch, 5 * zHalf, phi1, phi2); // a bit longer

      // the punched subtraction solids can be asymmetric and therefore have to be created twice:
      // one time in the "right" way, another time in the "reverse" way, because the "mirroring"
      // rotation around the y-axis will not only exchange +z and -z, but also +x and -x
      SubtractionSolid finalSolid0( wholeSolid, punchSolid, punchTransformer);
      SubtractionSolid finalSolid1( wholeSolid, punchSolid, punchTransmirror);

      // tube consists of vacuum (will later have two different daughters)
      Volume tubeLog0( volName + "_0", finalSolid0, material );
      Volume tubeLog1( volName + "_1", finalSolid1, material );
      if (isSensitive) {
        tubeLog0.setSensitiveDetector(sens);
        tubeLog1.setSensitiveDetector(sens);
      }
      tubeLog0.setVisAttributes(description, xmlMask.visStr() );
      tubeLog1.setVisAttributes(description, xmlMask.visStr() );

      // placement of the tube in the world, both at +z and -z
      PlacedVolume placed0 = envelope.placeVolume( tubeLog0, placementTransformer );
      PlacedVolume placed1 = envelope.placeVolume( tubeLog1, placementTransmirror );

      placed0.addPhysVolID("side",  1);
      placed1.addPhysVolID("side", -1);
      placed0.addPhysVolID("layer", counter);
      placed1.addPhysVolID("layer", counter);

      break;
    }

    case ODH::kUpstreamTrapezoid: {
      int maskSide = 1;
      if (xmlSection.hasAttr(_Unicode(side))) {
          maskSide = xmlSection.attr<int>(_Unicode(side));
      }

      double epsilon = 0.0;
      double extraThick = 0.0;
      if (xmlSection.hasAttr(_Unicode(extra_thickness))) {
          extraThick = xmlSection.attr<double>(_Unicode(extra_thickness));
          epsilon = xmlSection.attr<double>(_Unicode(epsilon));
      }

      // Trapezoid variables from xml
      const double origBaseHalf   = zHalf;         
      const double maskWidthHalf  = rInnerStart;   // rMin1: SeparatedBeamPipe_rmax
      const double offset         = rInnerEnd;     // rMin2: Distance from beam axis
      const double maskTopHalf    = rOuterStart;   // rMax1: Flat inner section half-length
      const double origHeightHalf = rOuterEnd;     // rMax2: Original thickness half-length

      // Slope-preserving math for the oversized trapezoid
      // As height increases by extraThick, the base expands proportionally so the slope angle remains invariant.
      const double maskHeightHalf = origHeightHalf + (extraThick / 2.0);
      const double maskBaseHalf   = maskTopHalf + (origBaseHalf - maskTopHalf) * (maskHeightHalf / origHeightHalf);

      Trd2 oversizedMask(maskBaseHalf, maskTopHalf, maskWidthHalf, maskWidthHalf, maskHeightHalf);

      // Determine Orientation and Position in the branch local frame
      Rotation3D maskOrientation;
      Position maskLocalPos;

      if (maskSide == 1) {
          // RIGHT MASK: Flat top facing inwards (-X)
          // Rotation in -90 degrees to make the trapezoid to face inwards
          maskOrientation = Rotation3D(RotationY(-90.0 * units::degree));
          maskLocalPos = Position(offset + maskHeightHalf, 0.0, 0.0);
      } else {
          // LEFT MASK: Flat top facing inwards (+X)
          // Rotation in 180 degrees to "mirror" the right trapezoid 

          maskOrientation = Rotation3D(RotationZ(180.0 * units::degree) * RotationY(-90.0 * units::degree));
          maskLocalPos = Position(-(offset + maskHeightHalf), 0.0, 0.0);
      }
      
      // Create the Transformation mapping from the mask's frame to the branch's frame
      Transform3D maskTransform(maskOrientation, maskLocalPos);

      // Boolean Subtraction to carve out the cylindrical wall
      // Multiply zHalf by 2.0 to guarantee the cutter is safely longer than the mask)
      // Multiplying outer wall to ensurre we subtracted all remaining trapezoid parts
      Tube outerTube(rInnerStart - epsilon, rInnerStart + 5 * maskHeightHalf, zHalf * 2.0);
      
      // To apply the subtraction, outerTube (which is at 0,0,0 in the branch frame) 
      // must be projected backward into the mask's local coordinate frame.
      Transform3D tubeInMaskFrame = maskTransform.Inverse();
      
      SubtractionSolid finalMaskSolid(oversizedMask, outerTube, tubeInMaskFrame);

      // Create the Logical Volume 
      Volume maskVol(volName, finalMaskSolid, material);
      if (isSensitive) {
          maskVol.setSensitiveDetector(sens);
      }
      maskVol.setVisAttributes(description, xmlMask.visStr());

      // Global Placements 
      Transform3D transformer, transmirror;
      if (rotationX) {
          transformer = Transform3D(RotationX(rotateAngle), RotateX(Position(0, 0, zPosition), rotateAngle));
          transmirror = Transform3D(RotationX(mirrorAngle), RotateX(Position(0, 0, zPosition), mirrorAngle));
      } else {
          transformer = Transform3D(RotationY(rotateAngle), RotateY(Position(0, 0, zPosition), rotateAngle));
          transmirror = Transform3D(RotationY(mirrorAngle), RotateY(Position(0, 0, zPosition), mirrorAngle));
      }

      // Place either into the host detector's enclosing volume (if <nest_in> was
      // given) or, as before, into our own assembly envelope.
      auto placeSection = [&](const Transform3D& volInWorld) -> PlacedVolume {
        if ( !hostVol.isValid() )
          return envelope.placeVolume(maskVol, volInWorld);

        if ( isSensitive )
          dd4hep::except("MDI_Mask",
                         "+++ %s: section '%s' is sensitive and requests nest_in. The "
                         "VolumeManager would scan it starting from the host's DetElement, so "
                         "its volID path loses %s's own 'system' field and the cellIDs would "
                         "be wrong. Drop either 'sensitive' or 'nest_in'.",
                         name.c_str(), volName.c_str(), name.c_str());

        const Transform3D volInHost = hostToWorld.Inverse() * volInWorld;
        ODH::NestingResult r =
          ODH::findDeepestContainer(hostVol, volInHost, maskVol.solid().ptr());

        if ( !r.mother.isValid() ) {
          dd4hep::printout(dd4hep::WARNING, "MDI_Mask",
                           "+++ %s: no volume inside '%s' contains '%s' - falling back to a "
                           "placement in '%s'. Geant4 will flatten the assemblies and report "
                           "an overlap.", name.c_str(), hostVol.name(), volName.c_str(),
                           envelope.name());
          return envelope.placeVolume(maskVol, volInWorld);
        }
        if ( r.samplesOutside > 0 )
          dd4hep::printout(dd4hep::WARNING, "MDI_Mask",
                           "+++ %s: %d of %d material sample points of '%s' lie outside the "
                           "selected host volume '%s' [%s] - the nesting is geometrically "
                           "inconsistent.", name.c_str(), r.samplesOutside, r.samplesTested,
                           volName.c_str(), r.mother.name(), r.path.c_str());
        dd4hep::printout(dd4hep::INFO, "MDI_Mask",
                         "+++ %s: nesting '%s' into '%s' (depth %d, path %s)",
                         name.c_str(), volName.c_str(), r.mother.name(), r.depth,
                         r.path.c_str());
        return r.mother.placeVolume(maskVol, r.childInMother);
      };

      // Multiply the global branch matrix by the local mask placement matrix
      Transform3D transformFwd = transformer * maskTransform;
      PlacedVolume placedFwd = placeSection(transformFwd);
      
      Transform3D transformBwd = transmirror * maskTransform;
      PlacedVolume placedBwd = placeSection(transformBwd);
      
      // Physical IDs
      placedFwd.addPhysVolID("side", 1).addPhysVolID("layer", counter).addPhysVolID("module", maskSide);
      placedBwd.addPhysVolID("side", -1).addPhysVolID("layer", counter).addPhysVolID("module", maskSide);

      break;
    }
    
    default: {
      throw std::runtime_error( " Mask_o1_v01_geo.cpp : fatal failure !! ??  " ) ;
    }

    }//end switch
    
  }//for all xmlSections

  //--------------------------------------
  Volume mother =  description.pickMotherVolume( tube ) ;
  PlacedVolume pv(mother.placeVolume(envelope));
  pv.addPhysVolID( "system", xmlMask.id() ) ; //.addPhysVolID("side", 0 ) ;

  tube.setVisAttributes( description, xmlMask.visStr(), envelope );

  tube.setPlacement(pv);

  return tube;
}
DECLARE_DETELEMENT(MDI_Mask_o1_v02,create_detector)
