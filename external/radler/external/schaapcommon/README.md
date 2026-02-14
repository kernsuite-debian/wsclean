# About schaapcommon
This repository contains utilities that are shared among different packages in the schaap-stack, including DP3, WSClean and EveryBeam.

## Requirements
When compiling `schaapcommon` as a stand-alone (static) library, the `aocommon` headers need to be available. `aocommon` can be cloned from https://gitlab.com/aroffringa/aocommon. To include the headers in the (cmake) build process, use the `AOCOMMON_INCLUDE_DIR` variable.

A `cmake` command typically reads

```
cmake -DAOCOMMON_INCLUDE_DIR=[PATH_TO_AOCOMMON/aocommon/include] -DCMAKE_INSTALL_PREFIX=[INSTALL_PATH] [PATH_TO_SCHAAPCOMMON]
```

# Dependants
The following repositories use schaapcommon:
- [DP3](https://git.astron.nl/RD/DP3)
- [EveryBeam](https://git.astron.nl/RD/EveryBeam)
- [Radler](https://git.astron.nl/RD/radler)
- [WSClean](https://gitlab.com/aroffringa/wsclean)
Note that these repositories thus depend on schaapcommon and it would be good practise to update them in the case breaking changes are made to this repository.
