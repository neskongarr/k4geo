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
#ifndef DDDETECTORS_ELLIPTOCONEHELPERS_H
#define DDDETECTORS_ELLIPTOCONEHELPERS_H 1

#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"

#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ECH { // EllipToConeHelpers

// Section types of EllipToCone_Mask_o1_v01 and EllipToCone_Beampipe_o1_v01. Both are elliptical cones centred on the
// z-axis (like ODH::kCenter), with semi-axes (a, b) along (x, y) changing linearly from the start to the end of the
// section. They differ only in how the polygonal cone is built, without any boolean solid:
typedef enum {
  kGeneralTraps = 0, // nSect general trapezoids (EightPointSolid, the G4GenericTrap of DD4hep), one per sector
  kTessSolids = 1,   // one tessellated solid
} EShapeType;

static EShapeType getShapeType(std::string const& type) {

  std::map<std::string, ECH::EShapeType> ShapeTypes;
  ShapeTypes["GeneralTraps"] = ECH::kGeneralTraps;
  ShapeTypes["TessSolids"] = ECH::kTessSolids;

  std::map<std::string, ECH::EShapeType>::const_iterator st = ShapeTypes.find(type);
  if (st == ShapeTypes.end()) {
    throw std::runtime_error("Unknown Shape Type for this geometry");
  }
  return st->second;
}

static bool checkForSensibleGeometry(double crossingAngle, EShapeType shapeType) {
  // both shape types are centred on the z-axis, so any crossing angle is fine
  switch (shapeType) {
  case kGeneralTraps:
  case kTessSolids:
    return true;
  default:
    printout(dd4hep::ERROR, "Mask/Beampip", "Unknown shape type, crossing angle %f.", crossingAngle);
    return false; // premature exit, dd4hep will abort now
  }
}

static double getCurrentAngle(double /*crossingAngle*/, EShapeType shapeType) {
  double tmpAngle;
  switch (shapeType) {
  case kGeneralTraps:
  case kTessSolids:
  default:
    tmpAngle = 0; // on the z-axis
    break;
  }

  return tmpAngle;
}

// General trapezoids.
// The vertices are points of the ellipse at fixed polar angles phi_i = 2 pi i / nSect. Every radial face of a sector is
// then a plane through the z-axis, so a wall sector lies exactly inside the tube sector of the same index, which is
// what allows the tube to be the mother of the wall without any boolean solid.
// The price: an edge at fixed polar angle is not a straight line of the cone when a/b changes along z, so the sectors
// lie slightly inside the true cone in the middle of the section (the end faces are exact).

// Point of the ellipse with semi-axes (a, b) at polar angle phi.
inline std::pair<double, double> ellipsePoint(double a, double b, double phi) {
  const double c = std::cos(phi), s = std::sin(phi);
  const double r = 1. / std::sqrt(c * c / (a * a) + s * s / (b * b));
  return {r * c, r * s};
}

// The wall: nSect hollow sectors between the inner (aMin, bMin) and the outer (aMax, bMax) ellipses, at -zHalf (1) and
// +zHalf (2). Vertex layout as in G4GenericTrap: 4 vertices at -zHalf, then 4 at +zHalf, clockwise in the xy-plane.
inline std::vector<dd4hep::Solid> constructGeneralTrapsWall(const std::string& solidName, double zHalf, int nSect,
                                                            double aMin1, double bMin1, double aMax1, double bMax1,
                                                            double aMin2, double bMin2, double aMax2, double bMax2) {
  std::vector<dd4hep::Solid> sectors(nSect);
  for (int i = 0; i < nSect; ++i) {
    const double phiCur = (2. * M_PI / nSect) * i;
    const double phiNxt = (2. * M_PI / nSect) * (i + 1);
    auto curIn1 = ellipsePoint(aMin1, bMin1, phiCur), nxtIn1 = ellipsePoint(aMin1, bMin1, phiNxt);
    auto curOut1 = ellipsePoint(aMax1, bMax1, phiCur), nxtOut1 = ellipsePoint(aMax1, bMax1, phiNxt);
    auto curIn2 = ellipsePoint(aMin2, bMin2, phiCur), nxtIn2 = ellipsePoint(aMin2, bMin2, phiNxt);
    auto curOut2 = ellipsePoint(aMax2, bMax2, phiCur), nxtOut2 = ellipsePoint(aMax2, bMax2, phiNxt);

    // clockwise: next inner, next outer, current outer, current inner
    double v[16] = {nxtIn1.first, nxtIn1.second, nxtOut1.first, nxtOut1.second,
                    curOut1.first, curOut1.second, curIn1.first, curIn1.second,
                    nxtIn2.first, nxtIn2.second, nxtOut2.first, nxtOut2.second,
                    curOut2.first, curOut2.second, curIn2.first, curIn2.second};

    sectors[i] = dd4hep::EightPointSolid(solidName + "_" + std::to_string(i), zHalf, v);
  }
  return sectors;
}

// The tube (including vacuum and wall): nSect solid sectors of the outer ellipse (aMax, bMax). The vertex on the
// z-axis is repeated, which makes the end faces triangles.
inline std::vector<dd4hep::Solid> constructGeneralTrapsTube(const std::string& solidName, double zHalf, int nSect,
                                                            double aMax1, double bMax1, double aMax2, double bMax2) {
  std::vector<dd4hep::Solid> sectors(nSect);
  for (int i = 0; i < nSect; ++i) {
    const double phiCur = (2. * M_PI / nSect) * i;
    const double phiNxt = (2. * M_PI / nSect) * (i + 1);
    auto cur1 = ellipsePoint(aMax1, bMax1, phiCur), nxt1 = ellipsePoint(aMax1, bMax1, phiNxt);
    auto cur2 = ellipsePoint(aMax2, bMax2, phiCur), nxt2 = ellipsePoint(aMax2, bMax2, phiNxt);

    // clockwise: next outer, current outer, axis, axis (the outer face is the one of the wall sector)
    double v[16] = {nxt1.first, nxt1.second, cur1.first, cur1.second, 0., 0., 0., 0.,
                    nxt2.first, nxt2.second, cur2.first, cur2.second, 0., 0., 0., 0.};

    sectors[i] = dd4hep::EightPointSolid(solidName + "_" + std::to_string(i), zHalf, v);
  }
  return sectors;
}

// Tessellated solids.
// The whole cone is one solid, so there are no faces between sectors and the vertices can be the exact points
// (a cos t, b sin t) of the cone: the straight line between the same t at -zHalf and +zHalf lies on the cone. Each
// twisted quadrilateral face becomes two triangles, with the same diagonal for the wall and the tube, so their outer
// surfaces coincide.

// The wall: outer side, inner side (normals towards the axis) and two annular end faces.
inline dd4hep::Solid constructTessSolidWall(const std::string& solidName, double zHalf, int nSect, double aMin1,
                                            double bMin1, double aMax1, double bMax1, double aMin2, double bMin2,
                                            double aMax2, double bMax2) {
  using Vertex = dd4hep::TessellatedSolid::Vertex;

  // rings of nSect vertices: outer at -zHalf, outer at +zHalf, inner at -zHalf, inner at +zHalf
  std::vector<Vertex> vertices;
  const double semiAxes[4][2] = {{aMax1, bMax1}, {aMax2, bMax2}, {aMin1, bMin1}, {aMin2, bMin2}};
  const double zRing[4] = {-zHalf, +zHalf, -zHalf, +zHalf};
  for (int ring = 0; ring < 4; ++ring) {
    for (int i = 0; i < nSect; ++i) {
      const double t = (2. * M_PI / nSect) * i;
      vertices.push_back(Vertex(std::cos(t) * semiAxes[ring][0], std::sin(t) * semiAxes[ring][1], zRing[ring]));
    }
  }
  const int o1 = 0, o2 = nSect, i1 = 2 * nSect, i2 = 3 * nSect;

  dd4hep::TessellatedSolid wall(solidName, vertices);
  for (int i = 0; i < nSect; ++i) {
    const int next = (i + 1) % nSect;
    wall.addFacet(o1 + i, o1 + next, o2 + next); // outer side
    wall.addFacet(o1 + i, o2 + next, o2 + i);
    wall.addFacet(i1 + i, i2 + next, i1 + next); // inner side
    wall.addFacet(i1 + i, i2 + i, i2 + next);
    wall.addFacet(o2 + i, o2 + next, i2 + next); // end face at +zHalf
    wall.addFacet(o2 + i, i2 + next, i2 + i);
    wall.addFacet(o1 + i, i1 + next, o1 + next); // end face at -zHalf
    wall.addFacet(o1 + i, i1 + i, i1 + next);
  }
  wall->CloseShape(true, true, false);
  return wall;
}

// The tube (including vacuum and wall): the solid cone of the outer ellipse (aMax, bMax).
inline dd4hep::Solid constructTessSolidTube(const std::string& solidName, double zHalf, int nSect, double aMax1,
                                            double bMax1, double aMax2, double bMax2) {
  using Vertex = dd4hep::TessellatedSolid::Vertex;

  // ring at -zHalf, ring at +zHalf, then the centres of the two end faces
  std::vector<Vertex> vertices;
  for (int i = 0; i < nSect; ++i) {
    const double t = (2. * M_PI / nSect) * i;
    vertices.push_back(Vertex(std::cos(t) * aMax1, std::sin(t) * bMax1, -zHalf));
  }
  for (int i = 0; i < nSect; ++i) {
    const double t = (2. * M_PI / nSect) * i;
    vertices.push_back(Vertex(std::cos(t) * aMax2, std::sin(t) * bMax2, +zHalf));
  }
  const int centre1 = vertices.size();
  vertices.push_back(Vertex(0., 0., -zHalf));
  const int centre2 = vertices.size();
  vertices.push_back(Vertex(0., 0., +zHalf));

  dd4hep::TessellatedSolid tube(solidName, vertices);
  for (int i = 0; i < nSect; ++i) {
    const int next = (i + 1) % nSect;
    tube.addFacet(i, next, nSect + next); // side, same triangles as the outer side of the wall
    tube.addFacet(i, nSect + next, nSect + i);
    tube.addFacet(centre1, next, i);                 // end face at -zHalf
    tube.addFacet(centre2, nSect + i, nSect + next); // end face at +zHalf
  }
  tube->CloseShape(true, true, false);
  return tube;
}

} // namespace ECH

#endif // DDDETECTORS_ELLIPTOCONEHELPERS_H
