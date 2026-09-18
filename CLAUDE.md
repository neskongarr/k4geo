# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

k4geo implements detector geometries in [DD4hep](https://dd4hep.web.cern.ch/) XML "compact files" plus C++
construction plugins, for use in the [Key4hep](https://github.com/key4hep) software stack. It ships two DD4hep
plugin libraries (`libk4geo.so` for geometry, `libk4geoG4.so` for Geant4-specific sensitive-detector/output
plugins) that are loaded by DD4hep at runtime based on the `type` attribute of `<detector>`/`<field>` elements
in compact XML.

Top-level layout is by experiment/collaboration: `FCCee/`, `FCChh/`, `CLIC/`, `ILD/`, `SiD/`, `MuColl/`,
`CaloTB/`, `FCalTB/`. Each contains a `compact/` tree of versioned detector XML (e.g.
`FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml`). C++ construction code lives centrally under
`detector/{tracker,calorimeter,fcal,muonSystem,other,PID,CaloTB}/`, not inside the experiment directories.

## Build

Requires DD4hep (>=1.31), Geant4, ROOT — normally obtained via the Key4hep stack:
```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh   # nightly, or sw.hsf.org/key4hep/setup.sh for stable
```

```bash
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j4
```

To make the rest of the key4hep stack pick up this local checkout instead of the central CVMFS `k4geo`,
run `k4_local_repo` from the repo's top-level directory (sets `$K4GEO` to point here). Compact files use
`${K4GEO}` to resolve shared resources such as field maps under `fieldmaps/`.

Useful CMake options: `BUILD_TESTING` (default ON), `INSTALL_COMPACT_FILES` (copy XML compact files into the
install area), `INSTALL_BEAMPIPE_STL_FILES` (download CAD STL files for the detailed MDI beampipe).

## Running/visualizing a detector

```bash
ddsim -G -N 10 path/to/DETECTOR.xml
```
`-G` fires the default particle gun; `-N` sets the event count. Detector-specific steering lives in external
repos (FCC-config, CLDConfig, etc.), not in this repo — see README.md for the mapping per detector.

Compact files are plain DD4hep XML assembled via `<include ref="..."/>`. A "standalone" compact file for a
subsystem (e.g. `FCCee/MDI/compact/*/*_standalone.xml`) is a minimal top-level file for debugging that one
subsystem in isolation, pulling in a `_collection.xml` that includes the actual geometry/materials/field-map
pieces. Field maps referenced from XML (`<field type="FieldBrBz" filename="...">`) point at `.root` files
under `fieldmaps/`; a wrong or stale `${K4GEO}` (or a renamed/missing map file) is a common source of a DD4hep
lookup error when parsing such an include — check that the referenced filename actually exists under
`$K4GEO/fieldmaps/` before assuming the XML structure is wrong.

## Tests

Tests are CTest-driven, defined in `test/CMakeLists.txt`, and mostly run `ddsim` against a compact file with a
particle gun (via the `add_ddsim_gun_test` helper) or check overlaps/name consistency. Run after building:
```bash
cd build
ctest                        # all tests
ctest -R IDEA_o1_v03          # tests matching a name
ctest -R IDEA_o1_v03 -V       # verbose output for one test
```
Tests set `k4geo_DIR` to the source tree (not an installed copy) so compact files resolve field maps from this
checkout, and fail on any `Exception`/`EXCEPTION`/`ERROR`/`Error` string in test output
(`FAIL_REGULAR_EXPRESSION` in `test/CMakeLists.txt`), so warnings that merely print "Error" in unrelated
context can fail a test — filter scripts (`test/scripts/filter_G4Tessellated_warnings.py`) exist for that.

Other test categories:
- `t_test_files_versions_FCCee`: checks that identically-named files under `FCCee/` are actually identical
  (exceptions listed in `utils/IdenticalFiles_ToBeIgnored.txt`).
- `t_name_check_*`: verifies the top-level `<detector>` name in a compact file matches its filename stem
  (`utils/detector_name_check.py`).
- `TestSensThickness`/`BeamCalZtest` (`test/src/`): standalone executables linked against `lcgeo` for
  geometric sanity checks on sensitive layer thickness.

## Adding/modifying detector construction code

- Construction plugins are registered with `DECLARE_DETELEMENT(<TypeName>, <create_function>)` (see e.g.
  `detector/tracker/DriftChamber_o1_v01.cpp`); `<TypeName>` is what compact XML's `<detector type="...">`
  must match to invoke that constructor.
- New `.cpp` sources must be added to the `file(GLOB sources ...)` list in the top-level `CMakeLists.txt`
  (per subdirectory), or they won't be compiled in.
- `detectorCommon/` and `detectorSegmentations/` are separate libraries with shared helpers/segmentation
  classes used across many detector constructors — check there before duplicating geometry utility code.
- `DCH_INFO_H_EXIST` gates drift-chamber code that depends on a DD4hep header only present from DD4hep 1.29+;
  code and tests behind it are conditionally compiled/registered.

## Conventions

- Detector version directories/files follow `<Detector>_o<N>_v<MM>` naming; CI checks (`detector_name_check.py`,
  the identical-files test) enforce naming and cross-copy consistency for FCCee in particular.
- `CODEOWNERS` assigns per-subsystem reviewers (e.g. `/FCCee/IDEA` → `@lopezzot`, `/FCCee/ALLEGRO` →
  `@BrieucF` et al.) — useful for figuring out who to ask about a given subdetector's intent.
- Formatting is enforced via pre-commit: `clang-format` for C++ (config in `.clang-format`), `ruff format`
  for Python (`.ruff.toml`), and XML well-formedness checking on all XML files.
