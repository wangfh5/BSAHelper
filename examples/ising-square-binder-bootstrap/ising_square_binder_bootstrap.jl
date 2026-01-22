using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using DataFrames
using JLD2
using LaTeXStrings
using Printf

using BSAHelper

const DATA_FILE = joinpath(@__DIR__, "data", "Ising-square-Binder.dat")
const RESULTS_DIR = joinpath(@__DIR__, "results")
const FIGS_DIR = joinpath(@__DIR__, "figs")

# Tuning parameter: the 2nd column of this dataset is 1/T (inverse temperature).
const X_COL = :beta
const Y_COL = :Binder
const Y_ERR_COL = :Binder_err

const CRITICAL_PARAM_NAME = "betac"
const ETA_TYPE = :none

const BOOTSTRAP_N = 100
const BOOTSTRAP_SEED = 42
const USE_CACHED_JLD = false

function read_ising_binder(path::AbstractString)
    isfile(path) || error("Data file not found: $path")

    L = Int[]
    beta = Float64[]
    binder = Float64[]
    binder_err = Float64[]

    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue

        cols = split(s)
        length(cols) < 4 && continue

        push!(L, parse(Int, cols[1]))
        push!(beta, parse(Float64, cols[2]))
        push!(binder, parse(Float64, cols[3]))
        push!(binder_err, parse(Float64, cols[4]))
    end

    df = DataFrame(L = L, beta = beta, Binder = binder, Binder_err = binder_err)
    return sort!(df, [:L, :beta])
end

mkpath(RESULTS_DIR)
mkpath(FIGS_DIR)

data = read_ising_binder(DATA_FILE)
L_values = sort(unique(data.L))

@printf("Loaded %d rows; L = %s\n", nrow(data), string(L_values))
@printf("β (= 1/T) range: %.6f - %.6f\n", minimum(data.beta), maximum(data.beta))

problem = BSAProblem(
    name = "2D Ising Binder ratio (square lattice)",
    L_values = L_values,
    data = data,
    x_col = X_COL,
    y_col = Y_COL,
    y_err_col = Y_ERR_COL,
    x_err_col = nothing
)

function find_bsa_binary()
    candidates = String[]
    env_bin = get(ENV, "BSA_BIN", "")
    isempty(env_bin) || push!(candidates, env_bin)
    push!(candidates, "bsa")
    push!(candidates, joinpath(homedir(), "bin", "bsa"))

    for candidate in candidates
        try
            return resolve_bsa_binary(candidate)
        catch
        end
    end

    return resolve_bsa_binary("bsa")
end

# BSA binary: try ENV["BSA_BIN"] → "bsa" in PATH → "$HOME/bin/bsa"
result_path = joinpath(RESULTS_DIR, "ising_square_binder_bootstrap.jld2")

bootstrap_result = nothing
plot_metadata = Dict{String,Any}()
plot_data_sections = Matrix{Float64}[]
phys_fmt = Dict{String,NamedTuple}()

if USE_CACHED_JLD && isfile(result_path)
    println("Found cached result; loading: $result_path")
    @load result_path bootstrap_result plot_metadata plot_data_sections phys_fmt
else
    bsa_binary = find_bsa_binary()

    # For bootstrap, disable MC (uncertainties come from bootstrap)
    bsa_cfg = BSAConfig(
        binary = bsa_binary,
        scaling_form = 0,
        use_mc = false,
        xscale = 1.0
    )

    # Binder ratio is dimensionless: fix c2 = 0, and fit only βc and ν (c1 = 1/ν)
    base_params = BSAParameters(
        Tc_init = 0.47,
        Tc_fixed = false,
        c1_init = 1.0,
        c1_fixed = false,
        c2_init = 0.0,
        c2_fixed = true,
        theta0_init = 1.0,
        theta1_init = 1.0,
        theta2_init = 1.0,
        theta_fixed = false
    )

    boot_cfg = BootstrapConfig(
        n_samples = BOOTSTRAP_N,
        seed = BOOTSTRAP_SEED,
        jitter_params = Dict(
            :Tc_init => 0.0,
            :c1_init => 0.0,
        ),
        y_sample_relative_error = 0.001,
        verbose = true
    )

    bootstrap_result = bootstrap_bsa_analysis(problem, boot_cfg, bsa_cfg, base_params)
    bootstrap_result === nothing && error("Bootstrap failed: bootstrap_bsa_analysis returned nothing")

    summary_path = joinpath(RESULTS_DIR, "bootstrap_summary.txt")
    save_bootstrap_summary(
        problem,
        boot_cfg,
        bootstrap_result,
        summary_path;
        critical_param_name = CRITICAL_PARAM_NAME,
        eta_type = ETA_TYPE
    )

    # Reconstruct the scaling function at bootstrap-mean parameters (for plotting and residuals)
    plot_cfg = BSAConfig(
        binary = bsa_binary,
        scaling_form = 0,
        use_mc = false,
        xscale = 1.8
    )

    plot_metadata, plot_data_sections, phys_fmt = prepare_bootstrap_plot_data(
        problem,
        bootstrap_result,
        plot_cfg;
        temp_dir = RESULTS_DIR,
        critical_param_name = CRITICAL_PARAM_NAME,
        eta_type = ETA_TYPE
    )

    # print chi2_red
    chi2red = plot_metadata["chi2"] / (plot_metadata["n_points"] - plot_metadata["n_freeparams"])
    println("chi2_red = $chi2red")

    @save result_path bootstrap_result plot_metadata plot_data_sections phys_fmt
    println("✓ Saved bootstrap result: $result_path")
end

isempty(plot_metadata) && error("Plot metadata is empty: prepare_bootstrap_plot_data failed")
isempty(phys_fmt) && error("Formatted physical quantities are empty: prepare_bootstrap_plot_data failed")

haskey(phys_fmt, "betac") && println("betac = $(phys_fmt["betac"].value_str) ± $(phys_fmt["betac"].error_str)")
haskey(phys_fmt, "nu") && println("ν  = $(phys_fmt["nu"].value_str) ± $(phys_fmt["nu"].error_str)")
haskey(phys_fmt, "c2") && println("c2 = $(phys_fmt["c2"].value_str) ± $(phys_fmt["c2"].error_str)")

## -------------------------------------------------------------------------- ##
##                        Plot the data collapse figure                       ##
## -------------------------------------------------------------------------- ##

using PyPlot
import BSAHelper: BSAPlotting

BSAPlotting.plot_bsa_data_collapse(
    plot_metadata,
    plot_data_sections,
    phys_fmt,
    FIGS_DIR;
    save_prefix = "ising_square_binder_bootstrap",
    critical_param_name = CRITICAL_PARAM_NAME,
    eta_type = ETA_TYPE,
    observable_label = L"U_4"
)
