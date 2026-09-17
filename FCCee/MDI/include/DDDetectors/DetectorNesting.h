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
//
// Geometric resolution of a mother volume inside another detector's volume tree.
//
// Some MDI components are physically located inside a volume that belongs to a
// different subdetector -- e.g. the synchrotron radiation mask sits in the bore of
// the beam pipe, i.e. inside the beam pipe's (solid) vacuum cone. Placing such a
// component into its own detector envelope makes it a *sibling* of that cone once
// Geant4 dissolves the DD4hep assemblies, and Geant4 then reports
//
//   "... apparently fully encapsulating volume ... at the same level!"
//
// DD4hep cannot express this relationship: declareParent() / <detector parent="...">
// / <composite> all resolve to the parent detector's own envelope volume (which is
// itself an assembly, so it does not help), and they allow only one mother per
// detector name -- while a crossing-angle geometry needs two, one per beam branch.
//
// This header resolves the mother *geometrically* instead, so that no internal
// volume name ever has to be written down: given the child's transform, it walks the
// host detector's placed-volume tree and returns the deepest real volume that
// contains the child, together with the child's transform in that volume's frame.
//
//==========================================================================
#ifndef DDDETECTORS_DETECTORNESTING_H
#define DDDETECTORS_DETECTORNESTING_H

#include "DD4hep/DetElement.h"
#include "DD4hep/DetectorTools.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Volumes.h"

#include <TGeoBBox.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoShape.h>
#include <TGeoVolume.h>

#include <algorithm>
#include <string>
#include <vector>

namespace ODH {

  //--------------------------------------------------------------------------
  // Matrix conversion
  //
  // NOTE: do *not* use dd4hep::detail::matrix::_transform(const TGeoMatrix*) for the
  // TGeo -> Transform3D direction. It is unsuitable here for two independent reasons
  // (DDCore/src/MatrixHelpers.cpp):
  //   1. it returns an all-zero (singular) rotation block whenever IsRotation() is
  //      false, which is exactly the case for the TGeoTranslation / TGeoIdentity
  //      nodes used for daughters placed at Position() or at identity;
  //   2. it scales translations by MM_2_CM (0.1 unless DD4HEP_USE_GEANT4_UNITS),
  //      whereas the forward path -- dd4hep::_addNode(..., const Transform3D&) in
  //      DDCore/src/Volumes.cpp -- writes pos.x()/y()/z() into the TGeoCombiTrans
  //      verbatim. The round trip would therefore shrink by a factor of ten.
  // The pair below is consistent with that forward path in both respects.
  //--------------------------------------------------------------------------

  /// Convert a Transform3D to a TGeoHMatrix, exactly as dd4hep::_addNode() does.
  inline TGeoHMatrix toTGeo(const dd4hep::Transform3D& tr) {
    dd4hep::Rotation3D rot;
    dd4hep::Position   pos;
    tr.GetRotation(rot);
    tr.GetTranslation(pos);
    double elements[9];
    rot.GetComponents(elements);
    TGeoRotation r;
    r.SetMatrix(elements);
    return TGeoHMatrix(TGeoCombiTrans(TGeoTranslation(pos.x(), pos.y(), pos.z()), r));
  }

  /// Convert a TGeoMatrix to a Transform3D. Unit- and identity-safe.
  /** GetRotationMatrix() is overridden by every concrete TGeo matrix class and
   *  returns the identity for pure translations, so no IsRotation() branch is needed.
   */
  inline dd4hep::Transform3D fromTGeo(const TGeoMatrix& m) {
    const Double_t* t = m.GetTranslation();
    const Double_t* r = m.GetRotationMatrix();
    return dd4hep::Transform3D(r[0], r[1], r[2], t[0],
                               r[3], r[4], r[5], t[1],
                               r[6], r[7], r[8], t[2]);
  }

  /// Transform mapping a detector's placement frame to the world frame.
  inline dd4hep::Transform3D placementToWorld(dd4hep::DetElement de) {
    dd4hep::detail::tools::PlacementPath path;
    dd4hep::detail::tools::placementPath(de, path);
    TGeoHMatrix m;
    dd4hep::detail::tools::placementTrafo(path, false, m); // local -> world
    return fromTGeo(m);
  }

  //--------------------------------------------------------------------------
  // Mother volume search
  //--------------------------------------------------------------------------

  /// Outcome of findDeepestContainer().
  struct NestingResult {
    /// Deepest real volume containing the child. Invalid if nothing contains it.
    dd4hep::Volume      mother;
    /// The child's transform expressed in `mother`'s frame.
    dd4hep::Transform3D childInMother;
    /// Node path from the search root to `mother`, for diagnostics.
    std::string         path;
    int                 depth          = 0;
    int                 samplesTested  = 0;
    int                 samplesOutside = 0;
  };

  namespace nesting_detail {

    struct Candidate {
      dd4hep::Volume vol;
      TGeoHMatrix    toRoot; ///< candidate frame -> search-root frame
      int            depth;
      std::string    path;
    };

    /// Collect every real volume below `current` whose shape contains `probeRoot`.
    /** `toRoot` maps `current`'s frame to the search-root frame; `probeRoot` is the
     *  probe point in the search-root frame.
     */
    inline void collect(dd4hep::Volume current, const TGeoHMatrix& toRoot,
                        const double probeRoot[3], int depth,
                        const std::string& path, std::vector<Candidate>& out) {
      TGeoVolume* vol = current.ptr();
      if (!vol)
        return;
      for (int i = 0, n = vol->GetNdaughters(); i < n; ++i) {
        TGeoNode*   node     = vol->GetNode(i);
        TGeoVolume* daughter = node ? node->GetVolume() : nullptr;
        if (!daughter)
          continue;

        TGeoHMatrix next(toRoot);
        next.Multiply(node->GetMatrix()); // daughter frame -> search-root frame
        const std::string sub = path + "/" + node->GetName();

        if (daughter->IsAssembly()) {
          // Assemblies are dissolved by the Geant4 converter, so one can never be the
          // answer -- but real volumes below it can. Descend unconditionally: an
          // assembly's Contains() only reports whether some daughter contains the
          // point, which the recursion establishes anyway.
          collect(daughter, next, probeRoot, depth + 1, sub, out);
          continue;
        }

        double probeLocal[3];
        next.MasterToLocal(probeRoot, probeLocal);
        if (!daughter->GetShape()->Contains(probeLocal))
          continue; // prunes the whole subtree

        out.push_back({dd4hep::Volume(daughter), next, depth + 1, sub});
        collect(daughter, next, probeRoot, depth + 1, sub, out);
      }
    }
  } // namespace nesting_detail

  /// Find the deepest real volume inside `root` that contains the given child.
  /** @param root           volume whose daughter tree is searched
   *  @param childInRoot    the child's transform in `root`'s frame
   *  @param childShape     optional: the child's solid, used for the sample check
   *  @param samplesPerAxis grid resolution of the sample check (<2 disables it)
   *
   *  The descent criterion is containment of the child's *origin*. A bounding-box
   *  test must not be used instead: TGeoSubtraction::ComputeBBox keeps the left
   *  operand's box, so a boolean child's box is the whole pre-subtraction solid and
   *  would falsely reject a perfectly valid mother. Candidates are therefore ranked
   *  deepest-first and validated with sample points taken from the real solid, so
   *  that a too-deep candidate is rejected in favour of a shallower one.
   */
  inline NestingResult findDeepestContainer(dd4hep::Volume             root,
                                            const dd4hep::Transform3D& childInRoot,
                                            const TGeoShape*           childShape     = nullptr,
                                            int                        samplesPerAxis = 5) {
    NestingResult     res;
    const TGeoHMatrix childMat = toTGeo(childInRoot);
    const double      probeRoot[3] = {childMat.GetTranslation()[0],
                                      childMat.GetTranslation()[1],
                                      childMat.GetTranslation()[2]};

    std::vector<nesting_detail::Candidate> candidates;
    nesting_detail::collect(root, TGeoHMatrix(), probeRoot, 0,
                            std::string(root.name()), candidates);
    if (candidates.empty())
      return res;

    std::stable_sort(candidates.begin(), candidates.end(),
                     [](const nesting_detail::Candidate& a,
                        const nesting_detail::Candidate& b) { return a.depth > b.depth; });

    // Fill `res` from a candidate and report whether the child's material fits in it.
    auto consider = [&](const nesting_detail::Candidate& c) {
      TGeoHMatrix inv = c.toRoot.Inverse();
      inv.Multiply(&childMat); // child frame -> candidate frame

      res.mother         = c.vol;
      res.depth          = c.depth;
      res.path           = c.path;
      res.childInMother  = fromTGeo(inv);
      res.samplesTested  = 0;
      res.samplesOutside = 0;

      const TGeoBBox* box = dynamic_cast<const TGeoBBox*>(childShape);
      if (!box || samplesPerAxis < 2)
        return true;

      const Double_t* o = box->GetOrigin();
      const double    d[3] = {box->GetDX(), box->GetDY(), box->GetDZ()};
      const int       N = samplesPerAxis;
      for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
          for (int iz = 0; iz < N; ++iz) {
            double p[3] = {o[0] + d[0] * (2.0 * ix / (N - 1) - 1.0),
                           o[1] + d[1] * (2.0 * iy / (N - 1) - 1.0),
                           o[2] + d[2] * (2.0 * iz / (N - 1) - 1.0)};
            if (!childShape->Contains(p))
              continue; // only real child material counts
            double q[3];
            inv.LocalToMaster(p, q);
            ++res.samplesTested;
            if (!c.vol.ptr()->GetShape()->Contains(q))
              ++res.samplesOutside;
          }
        }
      }
      return res.samplesOutside == 0;
    };

    for (const nesting_detail::Candidate& c : candidates) {
      if (consider(c))
        return res;
    }
    // Nothing fully contains the child: keep the deepest candidate and let the caller
    // report samplesOutside > 0.
    consider(candidates.front());
    return res;
  }

} // namespace ODH

#endif // DDDETECTORS_DETECTORNESTING_H
