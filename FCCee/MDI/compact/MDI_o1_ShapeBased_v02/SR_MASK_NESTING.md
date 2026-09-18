# Nesting the SR mask inside the beampipe vacuum

**Author:** nkontaxakis — 17/09/26

## Problem

`MDI_o1_ShapeBased_v02` places the Synchrotron Radiation mask (`SynchRadMask`, see
`Beampipe_o4_v06.xml`) *inside* the beampipe bore, whereas `MDI_o1_ShapeBased_v01` placed
it radially *outside* the wall (`rMin1 = SeparatedBeamPipe_rmax + 1mm + eps`), 
which is why v01 never hit this issue.

The `DD4hep_Beampipe_o1_v01` cpp driver 
builds each beampipe section as a **solid** vacuum cone spanning r = 0 → outer wall
radius, with the wall shell placed *inside* it as a daughter (`tube_UpStreamBeamPipe_1_0`
→ daughter `tube_UpStreamBeamPipe_1_wall_0`). 

The volumes defined in the `DD4hep_Mask_o1_v01` driver are however placed as
**siblings** to the vacuum cone, even though the mask geometrically sits *inside* the cone.
Running the repo's overlap macro (`k4geo/overlaps/overlap.mac` / `k4geo/utils/overlap.mac`) on the standalone MDI
geometry reported:

```
Overlap is detected for volume ...!BeamPipe_assembly_1#1!tube_UpStreamBeamPipe_1_0_8#8:8
apparently fully encapsulating volume ...!SynchRadMask_assembly_2#2!tube_SRmask_Right_0#0:0
at the same level!        ->  overlap at local point (7,-0.00039,-1289.81) by 9 mm
```

## Solution: geometric nesting via `<nest_in>`

Added a further key `<nest_in>` at detector level that resolves **geometrically** the mother volume 
of all detector `<sections>`, so no internal volume name is ever written in XML or C++:

```xml
<detector name="SynchRadMask" type="MDI_Mask_o1_v02" insideTrackingVolume="true" vis="TantalumVis">
    <parameter crossingangle="CrossingAngle" />
    <nest_in detector="BeamPipe"/>
    <!-- two <section type="UpstreamTrapezoid" side="±1" .../>, unchanged -->
</detector>
```

`<nest_in>` names a **detector** (`BeamPipe`), not a volume. At construction time the
driver:

1. Looks up the host `DetElement` by name (`Detector::detectors()`, a plain map —
   *not* `Detector::detector()`, whose path-search fallback throws an unhelpful
   "Unknown child" message for what is really an XML-ordering mistake). Requires the
   host to already be built (`<detector>` elements are processed in document order, and
   `BeamPipe` already precedes `SynchRadMask` in `Beampipe_o4_v06.xml`).
2. For each section placement, walks the host's placed-volume tree (`TGeoNode`s,
   descending through assemblies but never stopping on one, since Geant4 dissolves
   them) and finds the **deepest real volume whose shape contains the mask's
   placement origin**, mapped down through the accumulated node matrices.
3. Places the mask solid into that volume, with the relative transform derived from
   the actual matrices rather than any hand-supplied offset.
4. Falls back to the legacy behaviour (place into its own assembly) with a `WARNING`
   if nothing contains the child — this keeps the driver safe for detectors or
   sections that don't use `<nest_in>` at all, and safe for the case where nesting
   doesn't apply.

The two mask sections (`side="1"`/`side="-1"`) each produce two placements (the ±z
crossing-angle branches), and the resolver picks a **different** mother for each —
`tube_UpStreamBeamPipe_1_0` vs `_1` — automatically, with no XML naming either one.

Why the placement-origin containment test is safe here, and why a bounding-box test
would not be: the mask solid is a boolean subtraction (oversized trapezoid minus a
tube cutting away everything at r ≥ `rInnerStart - epsilon`). `TGeoSubtraction`'s
`ComputeBBox()` keeps the *left* operand's box — i.e. the bounding box of the
oversized, pre-subtraction shape — which extends well outside the beampipe radius even
though the actual (post-subtraction) solid does not. A bbox-based test would therefore
falsely reject the correct mother. The resolver instead tests the placement origin
(which is real mask material, safely inside the vacuum bore and outside the wall
shell) and, as a secondary check, samples grid points restricted to the true solid via
`TGeoShape::Contains()`.

A sensitive section combined with `<nest_in>` is explicitly rejected at construction
time: re-homing a sensitive volume into a foreign detector's volume tree would make
the `VolumeManager` scan it starting from the *host's* DetElement, corrupting the
volID/cellID encoding.

## Files changed

- **`FCCee/MDI/DDDetectors/include/DDDetectors/DetectorNesting.h`** (new) — the
  geometric mother-volume resolver (`ODH::findDeepestContainer`) plus a
  `Transform3D` ⟷ `TGeoHMatrix` conversion pair. Do **not** replace the latter with
  `dd4hep::detail::matrix::_transform(const TGeoMatrix*)`: it returns an all-zero
  (singular) rotation for any non-rotation node (exactly what the wall shell and the
  envelopes use) and scales translations by a unit factor that does not match the
  forward placement path (`dd4hep::_addNode`, `DDCore/src/Volumes.cpp`), which writes
  translations verbatim. The header ships its own pair to match the forward path
  exactly.
- **`FCCee/MDI/DDDetectors/src/Mask_o1_v02_geo.cpp`** (was `Mask_o1_v01_geo.cpp`,
  copied from the upstream DD4hep `Mask_o1_v01_geo.cpp` driver) — added the
  `<nest_in>` resolution and rewired the `UpstreamTrapezoid` case's two placements
  through it. The trapezoid geometry math itself (solid construction,
  `transformer`/`transmirror`) is untouched.
  Registered under a **new** plugin name, `MDI_Mask_o1_v02` (was
  `DD4hep_Mask_o1_v01`). This is required, not cosmetic:
  `libDDDetectors.so` (built from the DD4hep checkout) already registers
  `DD4hep_Mask_o1_v01`, and Gaudi's `PluginService` keeps the *first* registration for
  a given name silently — a second driver under the same name would be resolved by
  library load order, not by which one is actually wanted. Only `SynchRadMask` uses
  the new type; `BeamPipeShield`, `BeamPipeShield_noRot` and
  `HOMAbsorbers_PunchThrough` (also in `Beampipe_o4_v06.xml` / `HOMAbsorber.xml`) stay
  on the upstream `DD4hep_Mask_o1_v01`, since they need neither the trapezoid case nor
  nesting.
- **`FCCee/MDI/DDDetectors/include/DDDetectors/OtherDetectorHelpers.h`** — unchanged
  copy of the upstream header (already carries `kUpstreamTrapezoid` and the
  `"UpstreamTrapezoid"` string mapping); kept here only because the driver includes it
  locally.
- **`CMakeLists.txt`** (top-level) — added the single new source file to k4geo's
  `file(GLOB sources ...)` list and a matching `target_include_directories` entry.
  Deliberately **one explicit file**, not a glob over the whole
  `FCCee/MDI/DDDetectors/` folder: that folder used to be a full mirror of DD4hep's
  `DDDetectors/` (~30 drivers) and building all of it would re-register every
  `DD4hep_*` plugin name a second time, hitting the same silent first-registration
  problem. The folder has been trimmed to the 3 files actually needed.
- **`Beampipe_o4_v06.xml`** — `SynchRadMask` changed from
  `type="DD4hep_Mask_o1_v01"` to `type="MDI_Mask_o1_v02"`, plus the new
  `<nest_in detector="BeamPipe"/>` line. No other detector in this file changed.
  The upstream beampipe section itself was **not** re-split around the mask (an
  intermediate iteration of this fix had split `UpStreamBeamPipe_1` into two pieces
  around the mask position, removing ~30 mm of copper wall and vacuum coverage); the
  mask now nests inside the single, continuous `UpstreamClippedFront` section as
  originally intended.

An intermediate, now-superseded approach (moving the trapezoid case into the Beampipe
driver and hardcoding `TGeoManager::GetVolume("tube_UpStreamBeamPipe_1_0")` plus a
hand-computed `center_z`) was tried and discarded; the corresponding
`*_CHANGED.xml`/`*_CHANGED.cpp` scratch files have been removed.

## Why this stays entirely inside k4geo

The DD4hep checkout under `DD4hep/DDDetectors/` (the one that historically built these
plugins) was intentionally left untouched. `FCCee/MDI/DDDetectors/` here is a working
copy created specifically so MDI-specific driver changes live in k4geo, not in the
DD4hep fork.

## Out of scope

Running the overlap check on the full `MDI_o1_ShapeBased_v02_standalone.xml` also
shows 5 pre-existing overlaps between `BeamPipeShield` (`TaShield`, a Tungsten cone
shell, z ≈ 1195–2180 mm) and `HOMAbsorbers_PunchThrough`/`HOMAbsorber`
(z ≈ 1198–1299 mm), up to ~2 cm deep. These are **not** a containment problem — an
85×75 mm iron box and a tungsten cone genuinely occupying overlapping space — so
nesting cannot fix them (the resolver would correctly find no enclosing volume). They
are present identically in the `MDI_o1_ShapeBased_v01` baseline and are unrelated to
the SR mask; they need a dimensional or material decision, not a placement fix, and
are intentionally left untouched by this change.

## Verification

- **Build.** `MDI_Mask_o1_v02` registers exactly once in `libk4geo.components`; no
  "factory already found" warning; `DD4hep_Mask_o1_v01` unaffected.
- **TGeo overlap/extrusion check**
  (`DD4hep/DDCore/python/bin/checkOverlaps.py -t 0.001`): **0** illegal
  overlaps/extrusions on the standalone v02 geometry.
- **Construction log** shows the four expected nesting decisions, e.g.:
  ```
  MDI_Mask INFO +++ SynchRadMask: nesting 'tube_SRmask_Right' into 'tube_UpStreamBeamPipe_1_0' (depth 1, ...)
  ```
  (depth 1 confirms it lands directly on the vacuum cone, not inside the wall shell).
- **Geant4 overlap check** (`k4geo/utils/overlap.mac` via `ddsim --runType run`):
  overlap count dropped from **13 → 5** on the v02 standalone geometry; all four
  `tube_SRmask_{Right,Left}` placements report `OK!`; the remaining 5 are exactly the
  out-of-scope TaShield/HOMAbsorber pairs above, unchanged.
- **Pose invariance.** Compared the mask's global transform (translation + full
  rotation matrix, to the tree's node-matrix precision) between this version and a
  reference build using the unmodified upstream driver with the mask in its own
  assembly (no nesting): **bit-identical**. Nesting only changes the mask's mother
  volume, not its position in space.
- **Regression.** `MDI_o1_ShapeBased_v01_standalone.xml` (untouched driver/type) still
  reports the same 5 pre-existing overlaps as its own baseline — the no-`<nest_in>`
  code path is unaffected.
- **Smoke test.** 30 muon events through the forward region (`ddsim -G -N 30`, θ ∈
  [0.1°, 3°]) produced no navigation/geometry errors.
