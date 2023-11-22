# Changelog

## v1.2.2 - 2023-11-22

<!-- Release notes generated using configuration in .github/release.yml at main -->
### What's Changed

#### New Features

- Make DGE plugin more robust to changing or deleting subsets by @jfoster17 in https://github.com/gluesolutions/glue-genes/pull/58

**Full Changelog**: https://github.com/gluesolutions/glue-genes/compare/v1.2.1...v1.2.2

## v1.2.1 - 2023-11-17

<!-- Release notes generated using configuration in .github/release.yml at main -->
### What's Changed

- Require glue-qt > 0.3 to get region display
- Update documentation by @jfoster17 in https://github.com/gluesolutions/glue-genes/pull/56

**Full Changelog**: https://github.com/gluesolutions/glue-genes/compare/v1.2.0...v1.2.1

## v1.2.0 - 2023-11-15

<!-- Release notes generated using configuration in .github/release.yml at main -->
### What's Changed

#### New Features

- Support for loading and analyzing spatial transcriptomics data
- Support for display regions on images

#### Bug Fixes

- Fix ome_zarr reader to work for any number of supported dims by @jfoster17 in https://github.com/gluesolutions/glue-genes/pull/48

#### Other Changes

- Add tifffile file reader by @jfoster17 in https://github.com/gluesolutions/glue-genes/pull/49
- Allow DEG plugin to run versus rest (all other obs) by @jfoster17 in https://github.com/gluesolutions/glue-genes/pull/53
- Change differential gene expression plugin to be live-updating by @jfoster17 in https://github.com/gluesolutions/glue-genes/pull/54

**Full Changelog**: https://github.com/gluesolutions/glue-genes/compare/v1.1.0...v1.2.0

## v1.1.0 - 2023-06-14

<!-- Release notes generated using configuration in .github/release.yml at main -->
### What's Changed

- Support for spatial transcriptomics
- Basic Zarr and OME-Zarr readers
- Spaceranger data importer
- Optional OpenSlide data loader
- New MultiResolutionData object for tiled data
- Remove need for custom startup action
- New SyncComponent to designate components that are get up-to-date with subsets
- Uses glue-for-glue-genes rather than glue-core

#### Other Changes

- Update installation recommendation by @jfoster17 in https://github.com/gluesolutions/glue-genes/pull/23
- Update QComboBox to work with Qt6
- Store uns metadata and deal with named (DataFrame) obsm arrays
- Improve tests
- Make var names unique in anndata
- Make subset listener more robust
- Bug-fixes and naming code improvements

**Full Changelog**: https://github.com/gluesolutions/glue-genes/compare/v1.0.1...v1.1.0

## v1.0 - Unreleased

Initial Release of glue genes as an integrated package.
