# BSAHelper

[![文档](https://img.shields.io/badge/docs-latest-blue.svg)](https://wangfh5.github.io/BSAHelper/)

A small Julia helper package for interfacing with the external Bayesian Scaling Analysis (BSA) executable.

## Requirements

- You must build/install BSA yourself: `https://github.com/KenjiHarada/BSA.git`.
- Point `BSAHelper` to the executable via `ENV["BSA_BIN"]` or `BSAConfig(binary=...)`.
- For now, this package depends on [DataProcessforDQMC.jl](https://github.com/wangfh5/DataProcessforDQMC.jl.git) for specific statistics/formatting utilities.
- `DataProcessforDQMC.jl` is unregistered; you may need to add it to your project environment via `Pkg.add(url="https://github.com/wangfh5/DataProcessforDQMC.jl.git")` if needed.
