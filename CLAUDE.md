# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BSAHelper is a Julia package that provides a high-level interface to the external Bayesian Scaling Analysis (BSA) executable for finite-size scaling (FSS) analysis in statistical physics. The package wraps BSA command-line invocations, parses outputs, performs bootstrap uncertainty estimation, and generates publication-quality data collapse plots.

**Key requirement**: Users must build/install the BSA executable separately from https://github.com/KenjiHarada/BSA.git and configure `ENV["BSA_BIN"]` or pass `BSAConfig(binary=...)`.

## Development Commands

### Testing
```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run tests with multiple threads (for bootstrap parallelization)
julia --project --threads 4 -e 'using Pkg; Pkg.test()'
```

### Documentation
```bash
# Build documentation locally
julia --project=docs docs/make.jl

# Preview documentation with local server
DOCS_PREVIEW=true DOCS_PORT=8275 julia --project=docs docs/make.jl
# Then open http://localhost:8275
```

### Running Examples
```bash
# Run the Ising bootstrap example (requires BSA binary)
cd examples/ising-square-binder-bootstrap
julia --threads 4 ising_square_binder_bootstrap.jl
```

## Architecture

### Three-Module Design

The package is organized into three independent modules with clear separation of concerns:

1. **`bsa_core.jl` (BSACore)**: Low-level BSA interface
   - Wraps BSA executable invocation via `run_bsa_analysis()`
   - Parses BSA output files (`.op` format) into metadata and data sections
   - Handles parameter mapping between input/output conventions
   - **Critical**: BSA uses different parameter orderings for `scaling_form=0` vs `scaling_form=1`
     - Form 0: `p[0]=Tc, p[1]=c1, p[2]=c2, p[3-5]=theta`
     - Form 1: `p[0]=Tc, p[1]=c1, p[2]=c3, p[3]=c2, p[4-8]=theta`
   - Provides three-layer output architecture:
     - Raw metadata: `{"p0", "p1", "p2", ...}` (ambiguous, form-dependent)
     - Parameter dict: `{"Tc", "c1", "c2", "c3", ...}` (unambiguous)
     - Physical params: `{"Uc"/"Tc", "nu", "eta_phi"/"eta_psi", ...}` (user-friendly)

2. **`bsa_bootstrap.jl` (BSABootstrap)**: Statistical uncertainty estimation
   - Implements parallel bootstrap resampling with `@threads`
   - Generates synthetic datasets by resampling with Gaussian noise
   - Runs BSA on each bootstrap sample and collects parameter distributions
   - Uses robust statistics (IQR-based outlier filtering) to handle occasional pathological fits
   - Provides `prepare_bootstrap_plot_data()` to inject bootstrap errors into BSA output for plotting
   - **Design note**: Bootstrap samples are written with relative errors (`y_sample_relative_error`) to avoid treating each sample as exact

3. **`bsa_plotting.jl` (BSAPlotting)**: Visualization via package extension
   - Implemented as a package extension (`ext/BSAHelperPyPlotExt.jl`) to avoid hard PyPlot dependency
   - Only loads when user has `PyPlot` and `LaTeXStrings` installed
   - Generates data collapse plots with residuals, scaling functions, and error bars
   - Supports two plot modes:
     - `:full` - Complete analysis view (data + scaling function + residuals + χ²)
     - `:simple` - Publication-ready view (data only + fitting window indicator)
   - Handles X-axis errors (optional) and flexible observable labeling

### Data Flow for Bootstrap Analysis

```
Raw data (DataFrame)
  ↓ BSAProblem (filter by L_values)
  ↓ bootstrap_bsa_analysis() [parallel loop]
    ├─ Resample with Gaussian noise
    ├─ Write temporary BSA input file
    ├─ Run BSA with randomized initial parameters
    └─ Parse output → parameter dict
  ↓ BootstrapResult (param_means, param_stds, samples)
  ↓ extract_physical_params() → physical quantities
  ↓ format_physical_params() → formatted strings with significant digits
  ↓ prepare_bootstrap_plot_data() → inject bootstrap errors into metadata
  ↓ plot_bsa_data_collapse() → publication figure
```

### Parameter Naming Conventions

The package uses a three-layer naming system to handle BSA's ambiguous output format:

1. **BSA output layer** (`p[0]`, `p[1]`, ...): Raw indices from BSA, meaning depends on `scaling_form`
2. **Unambiguous parameter layer** (`Tc`, `c1`, `c2`, `c3`, `theta0-4`): Form-independent technical names
3. **Physical quantity layer** (`Uc`/`Tc`/`Jc`, `nu`, `eta_phi`/`eta_psi`, `omega`): User-friendly physics names

**Key functions**:
- `extract_parameter_dict(metadata)`: Layer 1 → Layer 2 (handles form-dependent mapping)
- `extract_physical_params(parameter_dict)`: Layer 2 → Layer 3 (applies physics transformations)
- `get_param_to_pidx_mapping(scaling_form)`: Reverse mapping for injecting bootstrap errors

## Important Implementation Details

### BSA Executable Invocation

- BSA is invoked via `Cmd()` with stdout/stderr redirected to files
- Command structure: `bsa [-c] [-f 1] [-w xscale] input.dat mask1 value1 mask2 value2 ...`
- Masks: `1` = fit parameter, `0` = fix parameter
- The package handles all command construction in `build_bsa_command()`

### Bootstrap Parallelization

- Uses `Base.Threads.@threads` for parallel bootstrap sampling
- Thread-safe via `ReentrantLock` for shared result collection
- Each thread writes to separate temporary files (named with thread ID)
- Progress bar updates are thread-safe via `ProgressMeter.next!()`

### Robust Statistics

Bootstrap fits can occasionally converge to pathological local minima (e.g., `Tc` far outside data range). To handle this:
- `compute_statistics()` applies IQR-based outlier filtering (Tukey fence with `k=10.0`)
- Removes extreme tails while preserving legitimate broad distributions
- Uses `DataProcessforDQMC.iqr_fence_filter()` for consistent filtering logic

### Plotting Extension System

The plotting module uses Julia 1.9+ package extensions:
- Core package has no PyPlot dependency
- Extension loads automatically when user installs PyPlot + LaTeXStrings
- Stub module in `src/bsa_plotting.jl` exports function names
- Actual implementation in `ext/bsa_plotting.jl`

### Error Propagation

- **Single BSA fit**: Uses Monte Carlo errors from BSA (stderr in output)
- **Bootstrap analysis**: Uses bootstrap standard deviation as error estimate
- **Error-of-error**: Calculated as `σ/√(2(n-1))` for significant digit determination
- **ν transformation**: `ν = 1/c1` with error `Δν = Δc1/c1²`

## Common Pitfalls

1. **BSA binary not found**: Set `ENV["BSA_BIN"]` or pass `binary=` to `BSAConfig`
2. **Form-dependent parameter mapping**: Always use `extract_parameter_dict()` instead of directly accessing `metadata["p2"]`
3. **Bootstrap instability**: Some datasets (like Ising Binder ratio) are sensitive to initial parameters and jitter ranges
4. **X-axis transformation**: For R-dependent FSS (correlation ratio), use `c1=0` so BSA uses X directly without `L^(1/ν)` scaling
5. **Y-axis scaling**: BSA handles `L^c2` scaling internally via fitting parameters, don't pre-scale data
6. **Threading**: Bootstrap benefits significantly from `julia --threads N`, but ensure thread count doesn't exceed CPU cores

## Dependencies

- **Registered**: DataFrames, Interpolations, Logging, Printf, ProgressMeter, Random, Statistics
- **Unregistered**: DataProcessforDQMC (from https://github.com/wangfh5/DataProcessforDQMC.jl)
- **Optional (for plotting)**: PyPlot, LaTeXStrings (loaded via package extension)

## Julia Version

This package targets Julia 1.11 and uses features like package extensions (requires Julia ≥1.9).

## Documentation

Documentation is written in Chinese (中文) and built with Documenter.jl. The docs are organized by module:
- `bsa_core.md`: BSA interface and parameter mapping
- `bsa_bootstrap.md`: Bootstrap analysis workflow
- `bsa_plotting.md`: Plotting API and customization
- `api.md`: Complete API reference

Documentation is automatically deployed to GitHub Pages on push to main branch.
