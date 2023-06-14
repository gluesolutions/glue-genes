# Changelog

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
