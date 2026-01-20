module BSABootstrap

using DataFrames
using Random
using Statistics: mean, std, quantile, median
using Printf
using Logging
using Base.Threads
using ProgressMeter
using DataProcessforDQMC: round_error, format_value_error, iqr_fence_filter

# Import BSACore from parent scope (must be included before this module)
import ..BSACore

# Import BSACore.extract_physical_params to extend it with BootstrapResult method
import ..BSACore: extract_physical_params

export BSAProblem, BootstrapConfig, BootstrapResult,
       bootstrap_bsa_analysis, save_bootstrap_summary,
       success_rate, prepare_bootstrap_plot_data,
       format_physical_params, extract_and_format_physical_params
# Note: extract_physical_params is imported and extended, not re-exported here

## -------------------------------------------------------------------------- ##
##                 Basic Structs defining a bootstrap analysis                ##
## -------------------------------------------------------------------------- ##

"""
    struct BSAProblem

Describe one observable together with its data columns.

# Design Philosophy
- Y scaling (L^c2) is handled by BSA fitting parameters (c2_init, c2_fixed), not pre-scaling
- For R-dependent FSS, use c1=0 so that X = R is used directly without transformation
- User should prepare clean data externally before defining BSAProblem, we only filter the data by L_values and remove missing entries.
"""
Base.@kwdef struct BSAProblem
    name::String
    L_values::Vector{Int}
    data::DataFrame
    x_col::Symbol
    y_col::Symbol
    y_err_col::Symbol
    x_err_col::Union{Symbol,Nothing} = nothing
end

"""
    struct BootstrapConfig

Control bootstrap sampling behaviour and parameter randomisation.

# Fields
- `n_samples`: Number of bootstrap iterations
- `seed`: Random seed for reproducibility
- `jitter_params`: Dict mapping parameter names to jitter radii (e.g., Dict(:Tc_init => 0.1, :c2_init => 0.05))
- `verbose`: Print progress every 100 iterations
- `keep_failed`: Keep failed fits in statistics (default: false)
- `tempdir`: Directory for temporary files (auto-created if nothing)
"""
Base.@kwdef struct BootstrapConfig
    n_samples::Int = 1000
    seed::Union{Int,Nothing} = nothing
    jitter_params::Dict{Symbol,Float64} = Dict{Symbol,Float64}()
    verbose::Bool = false
    keep_failed::Bool = false
    tempdir::Union{Nothing,String} = nothing
end

"""
    struct BootstrapResult

Hold statistics extracted from bootstrap iterations.

Stores unambiguous BSA parameters from `bsa_core.jl` (Tc, c1, c2, c3, ...).
Physical interpretation (Uc/Tc/Jc, ν, η_ψ/η_φ) is done downstream
using `extract_physical_params()`.
"""
struct BootstrapResult
    param_means::Dict{String,Float64}               # Unambiguous parameters: "Tc", "c1", "c2", "c3", ...
    param_stds::Dict{String,Float64}                # Standard deviations
    n_success::Int                                  # Successful fits
    n_trials::Int                                   # Total bootstrap samples
    bootstrap_samples::Vector{Dict{String,Float64}} # All samples for custom analysis
end

## -------------------------------------------------------------------------- ##
##                              helper functions                              ##
## -------------------------------------------------------------------------- ##

"""
    success_rate(result)

Return the fraction of successful bootstrap fits.
"""
success_rate(result::BootstrapResult) = result.n_trials == 0 ? 0.0 : result.n_success / result.n_trials

"""
    prepare_temp_dir(cfg)

Create (or reuse) the temporary directory used for intermediate BSA files.
"""
function prepare_temp_dir(cfg::BootstrapConfig)
    if cfg.tempdir === nothing
        return mktempdir()
    else
        isdir(cfg.tempdir) || mkpath(cfg.tempdir)
        return cfg.tempdir
    end
end

"""
    cleanup_temp_dir(cfg, path)

Remove the temporary directory if it was created automatically.
"""
function cleanup_temp_dir(cfg::BootstrapConfig, path::String)
    if cfg.tempdir === nothing
        isdir(path) && rm(path; recursive=true, force=true)
    end
end

"""
    compute_statistics(samples)

Compute robust mean and standard deviation for each parameter appearing in samples.

Why robust?
- Some bootstrap fits can converge to pathological local minima (e.g. Tc far outside data range).
- A tiny number of extreme outliers can destroy the plain mean/std and make the
  reconstructed collapse (which fixes parameters to the point estimate) look terrible,
  even when the bulk of samples is well-behaved.

Strategy:
- Drop non-finite values (NaN/Inf)
- Apply a very loose IQR-based outlier filter per-parameter (Tukey fence with a large
  multiplier) to remove extreme tails while preserving legitimate broad distributions.
"""
function compute_statistics(samples::Vector{Dict{String,Float64}})
    param_means = Dict{String,Float64}()
    param_stds = Dict{String,Float64}()
    all_keys = Set{String}()
    for sample in samples
        union!(all_keys, keys(sample))
    end
    for key in sort(collect(all_keys))
        values = [sample[key] for sample in samples if haskey(sample, key) && isfinite(sample[key])]
        isempty(values) && continue

        # NOTE: DataProcessforDQMC.iqr_fence_filter returns a NamedTuple:
        # (; keep, filtered, removed_min, removed_max)
        res = iqr_fence_filter(values; k=10.0, min_n=8)
        values_used = res.filtered
        param_means[key] = mean(values_used)
        param_stds[key] = std(values_used)
    end
    return param_means, param_stds
end

"""
    extract_physical_params(result::BootstrapResult; 
                           critical_param_name::String="Uc",
                           eta_type::Symbol=:none)

Convert bootstrap results to user-friendly physical quantities (batch processing version).

This is a wrapper around `BSACore.extract_physical_params()` that handles bootstrap
statistics by converting means and stds into (value, error) tuples.

# Input
`result.param_means` and `result.param_stds` contain unambiguous parameters {"Tc", "c1", "c2", ...}

# Returns
Dict{String, Tuple{Float64,Float64}} mapping physical quantity names to (mean, std)

# Design
- Constructs a parameter_dict with (mean, std) tuples from bootstrap statistics
- Delegates the actual conversion logic to `BSACore.extract_physical_params()`
- This keeps the physical transformation logic in one place (bsa_core.jl)
- Bootstrap module focuses on batch processing of statistical results
"""
function extract_physical_params(result::BootstrapResult; 
                                critical_param_name::String="Uc",
                                eta_type::Symbol=:none)
    # Construct parameter_dict from bootstrap statistics
    parameter_dict = Dict{String,Tuple{Float64,Float64}}()
    for key in keys(result.param_means)
        parameter_dict[key] = (result.param_means[key], result.param_stds[key])
    end
    
    # Delegate to bsa_core for the actual physical conversion
    return BSACore.extract_physical_params(parameter_dict; 
                                          critical_param_name=critical_param_name,
                                          eta_type=eta_type)
end

"""
    format_physical_params(phys_num, n_success; fmt=:decimal)

Format physical quantities with error bars based on error-of-std.

This is a pure formatting function that takes already-extracted physical
quantities and applies error-of-std based significant digit determination.

# Arguments
- `phys_num`: Dict{String, Tuple{Float64, Float64}} from `extract_physical_params`
- `n_success`: Number of successful bootstrap samples (for error-of-std calculation)
- `fmt`: Output format (:decimal or :scientific)

# Returns
Dict{String, NamedTuple} where each entry contains:
- `value`: Float64 mean value
- `error`: Float64 std (error bar)
- `err_of_err`: Float64 error of the error
- `digits`: Int significant digits for error
- `value_str`: String formatted value (e.g., "3.846")
- `error_str`: String formatted error (e.g., "0.012")
"""
function format_physical_params(phys_num::Dict{String,Tuple{Float64,Float64}},
                                n_success::Int;
                                fmt::Symbol=:decimal)
    phys_fmt = Dict{String,NamedTuple}()
    
    for (key, (val, err)) in phys_num
        if n_success > 1 && err > 0
            err_of_err = err / sqrt(2 * (n_success - 1))
            rounded_err, sigdigits = round_error(err, err_of_err)
            val_str, err_str = format_value_error(val, rounded_err, sigdigits; format=fmt)
        else
            err_of_err = 0.0
            sigdigits = -1
            val_str = @sprintf("%.6f", val)
            err_str = @sprintf("%.6f", err)
        end
        
        phys_fmt[key] = (
            value      = val,
            error      = err,
            err_of_err = err_of_err,
            digits     = sigdigits,
            value_str  = val_str,
            error_str  = err_str,
        )
    end
    
    return phys_fmt
end

"""
    extract_and_format_physical_params(result::BootstrapResult;
                                       critical_param_name="Uc",
                                       eta_type=:none, fmt=:decimal)

Convenience function combining `extract_physical_params` and `format_physical_params`.

This is the typical entry point for Bootstrap analyses: extracts physical
quantities from BootstrapResult and formats them with error-of-std based
significant digits.

# Data Flow
```
BootstrapResult
  → extract_physical_params(result; critical_param_name, eta_type)
  → phys_num
  → format_physical_params(phys_num, n_success; fmt)
  → phys_fmt
```
"""
function extract_and_format_physical_params(result::BootstrapResult;
                                            critical_param_name::String="Uc",
                                            eta_type::Symbol=:none,
                                            fmt::Symbol=:decimal)
    phys_num = extract_physical_params(result;
                                       critical_param_name=critical_param_name,
                                       eta_type=eta_type)
    return format_physical_params(phys_num, result.n_success; fmt=fmt)
end

## -------------------------------------------------------------------------- ##
##                          Core processing functions                         ##
## -------------------------------------------------------------------------- ##

"""
    filter_problem_data(problem)

Filter raw data by L_values and remove missing entries. 
"""
function filter_problem_data(problem::BSAProblem)
    mask = Vector{Bool}(undef, nrow(problem.data))
    idx = 1
    for row in eachrow(problem.data)
        # Only check L_values and missing data
        valid = row.L in problem.L_values
        if valid
            valid = !any(ismissing, (row.L, row[problem.x_col], row[problem.y_col], row[problem.y_err_col]))
            if valid && problem.x_err_col !== nothing
                valid &= !ismissing(row[problem.x_err_col])
            end
        end
        mask[idx] = valid
        idx += 1
    end
    filtered = problem.data[mask, :]
    return sort(filtered, [:L, problem.x_col])
end

"""
    write_bootstrap_dataset(problem, L_vals, x_sample, y_sample, y_err, path)

Write one bootstrap realisation to disk in the input format expected by BSA.

# Arguments
- `problem`: BSAProblem defining the FSS problem
- `L_vals`: Vector of system sizes for each data point
- `x_sample`: Bootstrap resampled X values
- `y_sample`: Bootstrap resampled Y values (no pre-scaling, BSA handles L^c2)
- `y_err`: Y error bars (already resampled or original)
- `path`: Output file path for BSA input

# Design
X values are written directly without transformation.
For R-dependent FSS, use c1=0 in BSAParameters so that BSA uses X directly.
Y and Y_err are written directly without scaling,
as BSA fitting parameters (c2_init, c2_fixed) control the L^c2 scaling.
"""
function write_bootstrap_dataset(problem::BSAProblem,
                                 L_vals::Vector{Int},
                                 x_sample::Vector{Float64},
                                 y_sample::Vector{Float64},
                                 y_err::Vector{Float64},
                                 path::String)
    open(path, "w") do io
        println(io, "# Bootstrap sample")
        println(io, "# L\tX\tObservable\tError")
        prev_L = nothing
        for idx in eachindex(L_vals)
            L_val = L_vals[idx]
            x_val = x_sample[idx]        # Write X directly without transformation
            y_val = y_sample[idx]        # No scaling - BSA handles L^c2
            err_val = abs(y_err[idx])    # No scaling - consistent with y_val
            if prev_L !== nothing && L_val != prev_L
                println(io)
            end
            prev_L = L_val
            @printf(io, "%6d\t% .9e\t% .9e\t% .9e\n", L_val, x_val, y_val, err_val)
        end
    end
end

"""
    randomise_parameters(params, cfg)

Apply random jitter to all configured parameters before a bootstrap iteration.
Each parameter in `cfg.jitter_params` is perturbed by a random value in the range ±radius.
"""
function randomise_parameters(params::BSACore.BSAParameters,
                              cfg::BootstrapConfig)
    if isempty(cfg.jitter_params)
        return params
    end
    
    fields = fieldnames(BSACore.BSAParameters)
    kwargs = Dict{Symbol,Float64}()
    
    for (fld, radius) in cfg.jitter_params
        if radius <= 0
            continue
        end
        
        if fld ∉ fields
            @warn "Unknown jitter field $(fld); skipping"
            continue
        end
        
        current_value = getfield(params, fld)
        jitter = (rand() - 0.5) * 2 * radius  # Uniform in [-radius, +radius]
        kwargs[fld] = current_value + jitter
    end
    
    if isempty(kwargs)
        return params
    end
    
    return BSACore.update_parameters(params; kwargs...)
end

## -------------------------------------------------------------------------- ##
##                                Main Functions                              ##
## -------------------------------------------------------------------------- ##

"""
    bootstrap_bsa_analysis(problem, cfg, bsa_cfg, base_params)

Run bootstrap resampling for the supplied problem definition and return statistics.

Returns BootstrapResult containing raw BSA parameters. Use `extract_physical_params()`
to convert to physical quantities (Uc/Tc, ν, η).
"""
function bootstrap_bsa_analysis(problem::BSAProblem,
                                cfg::BootstrapConfig,
                                bsa_cfg::BSACore.BSAConfig,
                                base_params)
    data = filter_problem_data(problem)
    if isempty(data)
        @error "No data points available after filtering"
        return nothing
    end

    base_params_struct = BSACore.ensure_parameters(base_params)

    if cfg.seed !== nothing
        Random.seed!(cfg.seed)
    end

    tmp_dir = prepare_temp_dir(cfg)
    samples = Dict{String,Float64}[]
    successes = 0

    println("="^60)
    println("Bootstrap BSA Analysis: $(problem.name)")
    println("="^60)
    println("Number of data points: $(nrow(data))")
    println("Number of bootstrap samples: $(cfg.n_samples)")
    println("Using $(nthreads()) thread(s) for parallel computation")
    println("L values: $(problem.L_values)")
    println("="^60)

    L_vals = Int.(data.L)
    x_vals = Float64.(data[:, problem.x_col])
    y_vals = Float64.(data[:, problem.y_col])
    y_err = Float64.(data[:, problem.y_err_col])
    x_err = problem.x_err_col === nothing ? zeros(length(x_vals)) : Float64.(data[:, problem.x_err_col])

    # Pre-allocate thread-safe storage for results
    samples_lock = ReentrantLock()
    successes_atomic = Atomic{Int}(0)
    
    # Create progress bar
    progress = Progress(cfg.n_samples, desc="Bootstrap sampling: ", 
                       barglyphs=BarGlyphs("[=> ]"), barlen=50)
    
    try
        # Parallel bootstrap loop with progress bar
        @threads for n in 1:cfg.n_samples
            # Each thread gets its own random samples (thread-safe)
            x_sample = x_vals .+ x_err .* randn(length(x_vals))
            y_sample = y_vals .+ y_err .* randn(length(y_vals))
            params = randomise_parameters(base_params_struct, cfg)

            # Use thread ID to avoid file conflicts
            tid = threadid()
            temp_fss = joinpath(tmp_dir, "bootstrap_$(n)_t$(tid).dat")
            temp_op = joinpath(tmp_dir, "bootstrap_$(n)_t$(tid).op")
            temp_log = joinpath(tmp_dir, "bootstrap_$(n)_t$(tid).log")

            try
                # Bootstrap realization should have zero error (each sample is treated as exact observation)
                zero_err = zeros(length(y_sample))
                write_bootstrap_dataset(problem, L_vals, x_sample, y_sample, zero_err, temp_fss)
                success = BSACore.run_bsa_analysis(bsa_cfg, params, temp_fss, temp_op, temp_log; silent=true)
                if success && isfile(temp_op)
                    metadata, _ = BSACore.parse_bsa_output(temp_op)
                    param_dict_full = BSACore.extract_parameter_dict(metadata)
                    # Extract only values (without errors) for bootstrap sampling
                    param_dict = Dict{String,Float64}(k => v[1] for (k, v) in param_dict_full)
                    if !isempty(param_dict)
                        lock(samples_lock) do
                            push!(samples, param_dict)
                        end
                        atomic_add!(successes_atomic, 1)
                    elseif cfg.keep_failed
                        lock(samples_lock) do
                            push!(samples, param_dict)
                        end
                    end
                elseif cfg.keep_failed
                    lock(samples_lock) do
                        push!(samples, Dict{String,Float64}())
                    end
                end
            finally
                for f in (temp_fss, temp_op, temp_log)
                    isfile(f) && rm(f; force=true)
                end
            end
            
            # Update progress bar (thread-safe)
            next!(progress)
        end
    finally
        cleanup_temp_dir(cfg, tmp_dir)
    end
    
    successes = successes_atomic[]

    if isempty(samples)
        @error "No successful bootstrap samples"
        return nothing
    end

    param_means, param_stds = compute_statistics(samples)

    result = BootstrapResult(param_means, param_stds, successes, cfg.n_samples, samples)

    println("="^60)
    println("Bootstrap results summary (raw parameters)")
    println("="^60)
    for key in sort(collect(keys(result.param_means)))
        @printf("  %s = %.6e ± %.6e\n", key, result.param_means[key], result.param_stds[key])
    end
    @printf("\nSuccess rate: %.1f%%\n", success_rate(result) * 100)
    println("Use extract_physical_params() to convert to physical quantities (Uc, ν, η)")
    println("="^60 * "\n")

    return result
end

## -------------------------------------------------------------------------- ##
##                          Functions After Bootstrap                         ##
## -------------------------------------------------------------------------- ##

"""
    save_bootstrap_summary(problem, config, result, output_file; 
                          critical_param_name="Uc", eta_type=:none)

Persist bootstrap statistics to a comprehensive human-readable summary.

Integrates context from `BSAProblem`, `BootstrapConfig`, and `BootstrapResult` 
to provide a self-contained summary including:
- Analysis metadata (problem name, L range, bootstrap config)
- Raw parameters (Tc, c1, c2, ...)
- Physical quantities (Uc/Tc/Jc, ν, η) with user-defined interpretation
- Success rate and distribution statistics

# Arguments
- `problem`: BSAProblem providing context (name, L_values)
- `config`: BootstrapConfig providing settings (n_samples, jitter_params)
- `result`: Bootstrap results containing param_means, param_stds, samples
- `output_file`: Path to save the summary
- `critical_param_name`: Physical name for critical point (default: "Uc")
- `eta_type`: Interpretation of c2 (:none, :eta_psi, :eta_phi)
"""
function save_bootstrap_summary(problem::BSAProblem, config::BootstrapConfig,
                                result::BootstrapResult, output_file::String;
                                critical_param_name::String="Uc",
                                eta_type::Symbol=:none)
    # Generate formatted physical quantities (contains value, error, value_str, error_str)
    phys_fmt = extract_and_format_physical_params(result;
                                                  critical_param_name=critical_param_name,
                                                  eta_type=eta_type)
    
    open(output_file, "w") do io
        println(io, "="^70)
        println(io, "Bootstrap FSS Analysis Summary")
        println(io, "="^70)
        println(io, "")
        
        # === Analysis Context ===
        println(io, "[Analysis Context]")
        println(io, "Observable: $(problem.name)")
        println(io, "L values: $(problem.L_values)")
        Lmin, Lmax = extrema(problem.L_values)
        println(io, "L range: $Lmin - $Lmax")
        println(io, "Bootstrap samples: $(config.n_samples)")
        println(io, "Success rate: $(round(success_rate(result) * 100, digits=1))% ($(result.n_success)/$(result.n_trials))")
        if !isempty(config.jitter_params)
            println(io, "Jitter parameters: $(config.jitter_params)")
        end
        println(io, "")
        
        # === Raw Parameters ===
        println(io, "[Parameters]")
        for key in sort(collect(keys(result.param_means)))
            mean_val = result.param_means[key]
            std_val = result.param_stds[key]
            @printf(io, "  %s = %.6e ± %.6e\n", key, mean_val, std_val)
        end
        println(io, "")
        
        # === Distribution Statistics ===
        println(io, "[Distribution Statistics]")
        @printf(io, "%-12s  %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s\n",
                "Parameter", "Mean(robust)", "Std(robust)", "Min", "25%", "Median", "75%", "Max", "Lower outliers", "Upper outliers")
        for key in sort(collect(keys(result.param_means)))
            values_raw = [sample[key] for sample in result.bootstrap_samples if haskey(sample, key) && isfinite(sample[key])]
            isempty(values_raw) && continue

            res = iqr_fence_filter(values_raw; k=10.0, min_n=8)
            values_used = res.filtered
            lower_outliers = res.removed_min
            upper_outliers = res.removed_max
            mean_val = mean(values_used)
            std_val = std(values_used)
            min_val = minimum(values_raw)
            q25 = quantile(values_raw, 0.25)
            median_val = median(values_raw)
            q75 = quantile(values_raw, 0.75)
            max_val = maximum(values_raw)
            @printf(io, "%-12s  %15.6e  %15.6e  %15.6e  %15.6e  %15.6e  %15.6e  %15.6e  %15d  %15d\n",
                    key, mean_val, std_val, min_val, q25, median_val, q75, max_val, lower_outliers, upper_outliers)
        end
        println(io, "")
        
        # === Physical Quantities ===
        if !isempty(phys_fmt)
            println(io, "[Physical Quantities]")
            for key in sort(collect(keys(phys_fmt)))
                entry = phys_fmt[key]
                @printf(io, "  %s = %s ± %s\n", key, entry.value_str, entry.error_str)
            end
            println(io, "")
        end
        println(io, "="^70)
        println(io, "Note: Physical quantities interpreted with critical_param_name=\"$critical_param_name\", eta_type=:$eta_type")
        println(io, "="^70)
    end
    println("✓ Bootstrap summary saved to: $output_file")
end

"""
    prepare_bootstrap_plot_data(problem, result, bsa_cfg; temp_dir=nothing,
                                critical_param_name="Uc", eta_type=:none)

Prepare Bootstrap analysis results for plotting with `plot_bsa_data_collapse`.

This function:
1. Runs BSA with all parameters fixed to bootstrap means
2. **Injects Bootstrap errors into metadata** (overwriting MC errors from single BSA fit)
3. **Generates formatted physical quantities** with proper significant digits based on error-of-std
4. Returns data ready for `plot_bsa_data_collapse`

# Workflow
1. Convert bootstrap mean parameters to `BSAParameters` (all fixed)
2. Run BSA with `bsa_cfg` (use_mc=false, parameters fixed)
3. Parse BSA output to get `(metadata, data_sections)`
4. **Overwrite parameters in metadata with Bootstrap mean ± std**
5. **Generate `phys_fmt`** via `extract_and_format_physical_params` for plotting

# Arguments
- `problem`: BSAProblem providing data and x_transform
- `result`: BootstrapResult containing param_means and param_stds (Bootstrap errors)
- `bsa_cfg`: BSAConfig for BSA settings (including xscale for plotting range)
- `temp_dir`: Temporary directory (auto-created if nothing)
- `critical_param_name`: Physical name for critical point (default: "Uc")
- `eta_type`: Interpretation of c2 (:none, :eta_psi, :eta_phi)

# Returns
Tuple of `(metadata, data_sections, phys_fmt)` where:
- `metadata::Dict`: Contains **Bootstrap errors** for parameters (Tc, c1, c2, ...)
- `data_sections::Vector{Matrix{Float64}}`:
  - `[1]`: Scaled data points
  - `[2]`: Scaling function [X, mu, sigma]
- `phys_fmt::Dict{String,NamedTuple}`: Formatted physical quantities with value_str/error_str

Returns `(Dict(), Vector{Matrix{Float64}}(), Dict{String,NamedTuple}())` on failure.

# Note
The returned `phys_fmt` should be passed to `plot_bsa_data_collapse` for proper
display of error bars with significant digits determined by error-of-std.
"""
function prepare_bootstrap_plot_data(
    problem::BSAProblem,
    result::BootstrapResult,
    bsa_cfg::BSACore.BSAConfig;
    temp_dir::Union{String,Nothing}=nothing,
    critical_param_name::String="Uc",
    eta_type::Symbol=:none
)
    # Prepare temporary directory
    cleanup_needed = (temp_dir === nothing)
    working_dir = cleanup_needed ? mktempdir() : temp_dir
    if !cleanup_needed && !isdir(working_dir)
        mkpath(working_dir)
    end
    
    try
        # Step 1: Extract bootstrap mean parameters (unambiguous: Tc, c1, c2, c3)
        param_means = result.param_means
        
        # Step 2: Build BSAParameters with all fitting params fixed to bootstrap means
        recon_params = BSACore.BSAParameters(
            Tc_init = param_means["Tc"],
            Tc_fixed = true,
            c1_init = param_means["c1"],
            c1_fixed = true,
            c2_init = get(param_means, "c2", 0.1),
            c2_fixed = haskey(param_means, "c2"),
            c3_init = get(param_means, "c3", 0.5),
            c3_fixed = haskey(param_means, "c3"),
            theta0_init = get(param_means, "theta0", 1.0),
            theta1_init = get(param_means, "theta1", 1.0),
            theta2_init = get(param_means, "theta2", 1.0),
            theta3_init = get(param_means, "theta3", 1.0),
            theta4_init = get(param_means, "theta4", 1.0),
            theta_fixed = true
        )
        
        # Step 3: Prepare data file (filter and transform)
        data = filter_problem_data(problem)
        if isempty(data)
            @error "No data points available after filtering"
            return (Dict(), Vector{Matrix{Float64}}(), Dict{String,NamedTuple}())
        end
        
        temp_fss = joinpath(working_dir, "bootstrap_recon.dat")
        L_vals = Int.(data.L)
        x_vals = Float64.(data[:, problem.x_col])
        y_vals = Float64.(data[:, problem.y_col])
        y_err = Float64.(data[:, problem.y_err_col])
        
        # Write data file (X values written directly without transformation)
        write_bootstrap_dataset(problem, L_vals, x_vals, y_vals, y_err, temp_fss)
        
        # Step 4: Run BSA with bsa_cfg (all parameters fixed, use_mc=false for reconstruction)
        temp_op = joinpath(working_dir, "bootstrap_recon.op")
        temp_log = joinpath(working_dir, "bootstrap_recon.log")
        
        recon_cfg = BSACore.BSAConfig(
            binary = bsa_cfg.binary,
            scaling_form = bsa_cfg.scaling_form,
            use_mc = false,  # No MC for reconstruction
            xscale = bsa_cfg.xscale  # Use xscale from bsa_cfg
        )
        
        success = BSACore.run_bsa_analysis(recon_cfg, recon_params, temp_fss, temp_op, temp_log; silent=false)
        
        if !success || !isfile(temp_op)
            @warn "Failed to reconstruct scaling function"
            return (Dict(), Vector{Matrix{Float64}}(), Dict{String,NamedTuple}())
        end
        
        # Step 5: Parse output and inject Bootstrap errors into metadata
        metadata, data_sections = BSACore.parse_bsa_output(temp_op)
        
        # Overwrite fitted parameters with Bootstrap mean ± std (key operation!)
        # This allows reusing plot_bsa_data_collapse with Bootstrap errors
        # 
        # CRITICAL: metadata uses "p0", "p1", "p2", ... as keys (from BSA output format)
        # We need to map param names ("Tc", "c1", "c2", ...) to correct p[i] indices
        scaling_form = get(metadata, "form", 0)
        param_to_pidx = BSACore.get_param_to_pidx_mapping(scaling_form)
        
        # Inject Bootstrap statistics into metadata using correct p[i] keys
        for (param_name, mean_val) in result.param_means
            std_val = get(result.param_stds, param_name, 0.0)
            if haskey(param_to_pidx, param_name)
                pidx_key = param_to_pidx[param_name]
                metadata[pidx_key] = (mean_val, std_val)  # Tuple format, same as parse_bsa_output
            end
        end
        
        # Step 6: Inject X errors to data_sections if available (BSA doesn't handle X errors, so we add them manually)
        if problem.x_err_col !== nothing && !isempty(data_sections)
            # Extract X errors from original problem data for each point in scaled_data
            scaled_data = data_sections[1]
            n_points = size(scaled_data, 1)
            x_errors = zeros(n_points)
            
            # Match each scaled data point back to original data by L and transformed X
            for i in 1:n_points
                X_scaled = scaled_data[i, 1]
                # L column depends on scaling_form:
                # - form=0: [X, Y, E, L, x, y, dy]
                # - form=1: [X1, Y, E, X2=1/L^c3, L, x, y, dy]
                L_val = scaling_form == 1 ? Int(scaled_data[i, 5]) : Int(scaled_data[i, 4])
                
                # Find matching row in original data
                matching_rows = data[(data.L .== L_val), :]
                if !isempty(matching_rows)
                    # Find closest X match (direct comparison, no transformation)
                    x_col_data = Float64.(matching_rows[:, problem.x_col])
                    idx = argmin(abs.(x_col_data .- X_scaled))
                    x_errors[i] = matching_rows[idx, problem.x_err_col]
                end
            end
            
            # Append X errors as the last column of scaled_data
            data_sections[1] = hcat(scaled_data, x_errors)
        end
        
        # Step 7: Generate formatted physical quantities with proper significant digits
        phys_fmt = extract_and_format_physical_params(result;
                                                      critical_param_name=critical_param_name,
                                                      eta_type=eta_type)
        
        return (metadata, data_sections, phys_fmt)
        
    finally
        # Cleanup if temp_dir was auto-created
        if cleanup_needed
            isdir(working_dir) && rm(working_dir; recursive=true, force=true)
        end
    end
end

end
