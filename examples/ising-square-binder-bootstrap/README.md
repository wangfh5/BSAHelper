# Minimal Example: 2D Ising (Square) Binder Ratio — Bootstrap + Plotting

This example demonstrates the full workflow of `BSAHelper` on a small, self-contained dataset:

1. Read raw Binder-ratio data from `data/Ising-square-Binder.dat`
2. Build a `BSAProblem`
3. Run `BSABootstrap.bootstrap_bsa_analysis` to estimate `Tc` and `ν`
4. Reconstruct the scaling function and plot a data-collapse figure via `BSAPlotting.plot_bsa_data_collapse`

## Dataset provenance

The file `data/Ising-square-Binder.dat` is copied from the official BSA repository sample:

- https://github.com/KenjiHarada/BSA (path: `CC2/Sample/Ising-square-Binder.dat`)

### Caveat: this bootstrap workflow may not be robust for all datasets

This example is intended as a *workflow demo*, not a guaranteed “drop-in” analysis pipeline for arbitrary datasets.
In particular, for the Ising-model Binder-ratio dataset used here, we have observed that the fitted results can be **not very stable** (sensitive to initial parameters and may be trapped in local minima).
The initial parameters, as well as their jittering ranges, has to be carefully chosen. 
We find that for this dataset, jittering the initial parameters is unstable. 

## Requirements

- Julia `1.11` (this package targets Julia 1.11)
- A working BSA executable (from the upstream repository)

### Build / install BSA

One typical way (macOS example, using the upstream `CC2` target):

```bash
git clone https://github.com/KenjiHarada/BSA.git
cd BSA/CC2
make
```

Then either:

- Point `BSAHelper` to the binary:

```bash
export BSA_BIN=/path/to/BSA/CC2/new_bfss
```

or:

- Put it on your PATH as `bsa`:

```bash
mkdir -p ~/bin
ln -s /path/to/BSA/CC2/new_bfss ~/bin/bsa
export PATH="$HOME/bin:$PATH"
```

## Get the code

This example is shipped with `BSAHelper`. If you haven't already, clone the repository and enter this folder:

```bash
git clone https://github.com/wangfh5/BSAHelper.git
cd BSAHelper/examples/ising-square-binder-bootstrap
```

## Run the example

From this directory:

```bash
julia ising_square_binder_bootstrap.jl
```

or use multiple threads to speed up the computation:

```bash
julia --threads 4 ising_square_binder_bootstrap.jl
```

The script activates and instantiates a local Julia environment (`Project.toml`) automatically on first run.
It will also fetch `DataProcessforDQMC.jl` (an unregistered dependency) as needed.
The first run may take longer because `PyPlot.jl`/`PyCall.jl` can set up a Python backend automatically.

## Outputs

After a successful run, you should see:

- `results/bootstrap_summary.txt`
- `results/ising_square_binder_bootstrap.jld2` (saved bootstrap samples + plot data)
- `figs/ising_square_binder_bootstrap_collapse.png`

## Re-run fast (without recomputing)

By default the script uses the cached JLD2 file if it exists.
To force recomputation, either:

- Delete `results/ising_square_binder_bootstrap.jld2`, or
- Set `USE_CACHED_JLD = false` inside `ising_square_binder_bootstrap.jl`

## Inspect bootstrap samples

You can load the saved JLD2 file and inspect individual bootstrap samples:

```julia
using JLD2
@load "results/ising_square_binder_bootstrap.jld2" bootstrap_result

bootstrap_result.bootstrap_samples  # Vector{Dict{String,Float64}}
```

## Expected results (order of magnitude)

For this 2D Ising (square lattice) dataset, the estimates should be close to:

- `Tc ≈ 0.4407` (this is `βc`)
- `ν ≈ 1.0`
