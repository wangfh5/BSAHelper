module BSACore

using Interpolations
using Printf

export BSAConfig, BSAParameters, run_bsa_analysis, resolve_bsa_binary
export ensure_parameters, update_parameters, get_param_to_pidx_mapping
export parse_bsa_output, extract_parameter_dict, extract_physical_params, print_summary
export chi2_interp, chi2red_interp

## -------------------------------------------------------------------------- ##
##                    Basic Structs defining a BSA calling                    ##
## -------------------------------------------------------------------------- ##

"""
    struct BSAConfig

Configuration for the external BSA executable (executable location and options).

`binary` can be:
- An absolute or relative path to the executable, or
- A command name available in PATH.

If `ENV["BSA_BIN"]` is set, it will be used as the default binary.
"""
Base.@kwdef struct BSAConfig
    binary::String = get(ENV, "BSA_BIN", "bsa")
    scaling_form::Int = 0
    use_mc::Bool = true
    xscale::Float64 = 1.0  # Range of scaling function output (1.0 = largest L range only)
end

"""
    struct BSAParameters

Parameter initial guesses and fixed flags passed to the BSA binary.
"""
Base.@kwdef struct BSAParameters
    Tc_init::Float64 = 3.85
    Tc_fixed::Bool = false
    c1_init::Float64 = 0.9
    c1_fixed::Bool = false
    c2_init::Float64 = 0.1
    c2_fixed::Bool = false
    c3_init::Float64 = 0.5
    c3_fixed::Bool = false
    theta0_init::Float64 = 1.0
    theta1_init::Float64 = 1.0
    theta2_init::Float64 = 1.0
    theta3_init::Float64 = 1.0
    theta4_init::Float64 = 1.0
    theta_fixed::Bool = false
end

const _BSA_PARAM_FIELDS = fieldnames(BSAParameters)

"""
    resolve_bsa_binary(binary::String) -> String

Resolve the BSA executable path.

Accepts either a filesystem path or a command in PATH. Throws an error with a
clear message if the binary cannot be found.
"""
function resolve_bsa_binary(binary::String)
    isempty(binary) && error("BSA binary path is empty. Set BSAConfig(binary=...) or ENV[\"BSA_BIN\"].")
    if isfile(binary)
        return binary
    end
    resolved = Sys.which(binary)
    resolved === nothing && error("BSA binary not found: $(binary). Set BSAConfig(binary=...) or ENV[\"BSA_BIN\"], or ensure it is in PATH.")
    return resolved
end

"""
    ensure_parameters(params)

Normalize parameter configuration to a `BSAParameters` instance.
"""
ensure_parameters(params::BSAParameters) = params
function ensure_parameters(params::NamedTuple)
    return BSAParameters(; params...)
end

"""
    update_parameters(params; kwargs...)

Return a copy of `params` with selected fields updated.
"""
function update_parameters(params::BSAParameters; kwargs...)
    data = Dict{Symbol,Any}()
    for field in _BSA_PARAM_FIELDS
        data[field] = getfield(params, field)
    end
    for (k, v) in kwargs
        data[k] = v
    end
    return BSAParameters(; data...)
end

"""
    get_param_to_pidx_mapping(scaling_form::Int)

Build reverse mapping from parameter names to BSA output indices.

This provides the inverse of the mapping used in `extract_parameter_dict()`,
allowing conversion from unambiguous parameter names back to p[i] keys.

# Parameter Mapping
- `scaling_form=0`: Tc→p0, c1→p1, c2→p2, theta0-2→p3-5
- `scaling_form=1`: Tc→p0, c1→p1, c3→p2, c2→p3, theta0-4→p4-8

# Returns
Dict{String,String} mapping parameter names to p[i] keys

# Example
```julia
mapping = get_param_to_pidx_mapping(0)
# {"Tc" => "p0", "c1" => "p1", "c2" => "p2", ...}
```

# Usage
Used primarily in `prepare_bootstrap_plot_data()` to inject Bootstrap statistics
into metadata using the correct p[i] keys that `extract_parameter_dict()` expects.
"""
function get_param_to_pidx_mapping(scaling_form::Int)
    param_to_pidx = Dict{String,String}()
    
    # Common mappings (same for both forms)
    param_to_pidx["Tc"] = "p0"
    param_to_pidx["c1"] = "p1"
    
    # Form-dependent mappings
    if scaling_form == 0
        param_to_pidx["c2"] = "p2"
        param_to_pidx["theta0"] = "p3"
        param_to_pidx["theta1"] = "p4"
        param_to_pidx["theta2"] = "p5"
    elseif scaling_form == 1
        param_to_pidx["c3"] = "p2"
        param_to_pidx["c2"] = "p3"
        param_to_pidx["theta0"] = "p4"
        param_to_pidx["theta1"] = "p5"
        param_to_pidx["theta2"] = "p6"
        param_to_pidx["theta3"] = "p7"
        param_to_pidx["theta4"] = "p8"
    else
        @error "Invalid scaling_form: $scaling_form (must be 0 or 1)"
    end
    
    return param_to_pidx
end

## -------------------------------------------------------------------------- ##
##                 Core functions for formatting a BSA command                ##
## -------------------------------------------------------------------------- ##

mask_flag(fixed::Bool) = fixed ? "0" : "1"

"""
    build_parameter_segment(cfg, params)

Construct the mask/value segment appended to the BSA command.
"""
function build_parameter_segment(cfg::BSAConfig, params::BSAParameters)
    theta_mask = mask_flag(params.theta_fixed)
    if cfg.scaling_form == 0
        return [
            mask_flag(params.Tc_fixed), string(params.Tc_init),
            mask_flag(params.c1_fixed), string(params.c1_init),
            mask_flag(params.c2_fixed), string(params.c2_init),
            theta_mask, string(params.theta0_init),
            theta_mask, string(params.theta1_init),
            theta_mask, string(params.theta2_init)
        ]
    elseif cfg.scaling_form == 1
        return [
            mask_flag(params.Tc_fixed), string(params.Tc_init),
            mask_flag(params.c1_fixed), string(params.c1_init),
            mask_flag(params.c3_fixed), string(params.c3_init),
            mask_flag(params.c2_fixed), string(params.c2_init),
            theta_mask, string(params.theta0_init),
            theta_mask, string(params.theta1_init),
            theta_mask, string(params.theta2_init),
            theta_mask, string(params.theta3_init),
            theta_mask, string(params.theta4_init)
        ]
    else
        error("Invalid scaling_form: $(cfg.scaling_form)")
    end
end

"""
    build_bsa_command(cfg, params, input_file)

Assemble the full command vector for invoking the BSA executable.
"""
function build_bsa_command(cfg::BSAConfig, params::BSAParameters, input_file::String)
    binary = resolve_bsa_binary(cfg.binary)
    cmd = String[binary]
    if cfg.use_mc
        push!(cmd, "-c")
    end
    if cfg.scaling_form == 1
        push!(cmd, "-f", "1")
    elseif cfg.scaling_form != 0
        error("Invalid scaling_form: $(cfg.scaling_form)")
    end
    push!(cmd, "-w", string(cfg.xscale))
    push!(cmd, input_file)
    append!(cmd, build_parameter_segment(cfg, params))
    return cmd
end

## -------------------------------------------------------------------------- ##
##                  Main functions for running a BSA analysis                 ##
## -------------------------------------------------------------------------- ##

"""
    run_bsa_analysis(cfg, params, input_file, output_op, output_log; silent=false)

Invoke the BSA executable with the provided configuration and parameter masks.
Returns true on success.
"""
function run_bsa_analysis(cfg::BSAConfig, params::BSAParameters,
                          input_file::String, output_op::String, output_log::String;
                          silent::Bool=false)
    if !isfile(input_file)
        @error "Input file not found: $input_file"
        return false
    end
    cmd_parts = build_bsa_command(cfg, params, input_file)
    if !silent
        println("Running BSA analysis...")
        println("Command: ", join(cmd_parts, " "))
    end
    try
        open(output_op, "w") do op_io
            open(output_log, "w") do log_io
                run(pipeline(Cmd(cmd_parts), stdout=op_io, stderr=log_io))
            end
        end
        if !silent
            println("✓ BSA analysis completed")
            println("  Output: $output_op")
            println("  Log: $output_log")
        end
        return true
    catch e
        if !silent
            @error "BSA analysis failed" exception=(e, catch_backtrace())
        end
        return false
    end
end

"""
    run_bsa_analysis(input_file, output_op, output_log; kwargs...)

Convenience overload accepting keyword arguments for configuration and parameters.
"""
function run_bsa_analysis(input_file::String, output_op::String, output_log::String;
                          kwargs...)
    cfg_keys = Set([:binary, :scaling_form, :use_mc, :xscale])
    param_keys = Set([:Tc_init, :Tc_fixed, :c1_init, :c1_fixed, :c2_init, :c2_fixed,
                      :c3_init, :c3_fixed, :theta0_init, :theta1_init, :theta2_init,
                      :theta3_init, :theta4_init, :theta_fixed])
    cfg_kwargs = Dict{Symbol,Any}()
    param_kwargs = Dict{Symbol,Any}()
    silent = get(Dict(kwargs), :silent, false)
    for (k, v) in kwargs
        if k == :silent
            continue
        elseif k in cfg_keys
            cfg_kwargs[k] = v
        elseif k in param_keys
            param_kwargs[k] = v
        else
            param_kwargs[k] = v
        end
    end
    cfg = BSAConfig(; cfg_kwargs...)
    params = BSAParameters(; param_kwargs...)
    return run_bsa_analysis(cfg, params, input_file, output_op, output_log; silent=silent)
end

## -------------------------------------------------------------------------- ##
##                 Functions for parsing a BSA analysis output               ##
## -------------------------------------------------------------------------- ##

"""
    parse_bsa_output(filename)

Parse a BSA .op file into metadata and numeric data sections.
"""
function parse_bsa_output(filename::String)
    lines = readlines(filename)
    metadata = Dict{String, Any}()
    data_sections = Vector{Matrix{Float64}}()
    current_section = Vector{Vector{Float64}}()
    blank_count = 0
    for line in lines
        if startswith(line, "# Number of data points")
            metadata["n_points"] = parse(Int, split(line, '=')[2])
        elseif startswith(line, "# Number of free parameters")
            metadata["n_freeparams"] = parse(Int, split(line, '=')[2])
        elseif startswith(line, "# Scaling form")
            metadata["form"] = parse(Int, split(split(line, ':')[2])[1])
        elseif startswith(line, "# chi^2 =")
            metadata["chi2"] = parse(Float64, split(line, '=')[2])
        elseif startswith(line, "# p[")
            idx = parse(Int, match(r"\[(\d+)\]", line).captures[1])
            parts = split(split(line, '=')[2])
            value = parse(Float64, parts[1])
            stderr = length(parts) > 1 ? parse(Float64, parts[2]) : NaN
            metadata["p$(idx)"] = (value, stderr)
        end
        if startswith(line, "#")
            continue
        end
        stripped = strip(line)
        if isempty(stripped)
            blank_count += 1
            if blank_count >= 2 && !isempty(current_section)
                push!(data_sections, hcat(current_section...)')
                current_section = Vector{Vector{Float64}}()
                blank_count = 0
            end
        else
            blank_count = 0
            try
                values = parse.(Float64, split(stripped))
                push!(current_section, values)
            catch
                continue
            end
        end
    end
    if !isempty(current_section)
        push!(data_sections, hcat(current_section...)')
    end
    return metadata, data_sections
end

"""
    chi2_interp(metadata, data_sections; sigma_xy=nothing) -> Float64

Recompute an approximate χ² by interpolating the scaling function and evaluating

    χ² = Σᵢ ((Yᵢ - F(Xᵢ)) / Eᵢ)²

This is a *naive* goodness-of-fit measure, meant to be consistent with the
residuals plot logic (linear interpolation with linear extrapolation).

# Inputs
- `data_sections[1]`: scaled data table (at least 3 columns: X, Y, E).
  If an extra column `xerr` is present (in the same convention as `BSAPlotting.plot_bsa_data_collapse`),
  this function will propagate X errors into an effective variance:
  `σ_eff^2 = σ_y^2 + (df/dx)^2 σ_x^2`.
- `data_sections[2]`: scaling function table (at least 2 columns: X_func, mu_func; sigma is ignored)

# Keywords
- `sigma_xy`: optional covariance between X and Y for each point (same coordinates as scaled_data),
  used only when `xerr` is present.
  When provided, the effective variance becomes:
  `σ_eff^2 = σ_y^2 + (df/dx)^2 σ_x^2 - 2 (df/dx) σ_xy`.
  - If `sigma_xy` is a Real, the covariance is assumed to be constant for all points.
"""
function chi2_interp(metadata::Dict{String,Any},
                     data_sections::Vector{<:AbstractMatrix};
                     sigma_xy::Union{AbstractVector,Real,Nothing}=nothing)
    length(data_sections) >= 2 || throw(ArgumentError("data_sections must contain at least 2 sections: scaled_data and scaling_func"))

    scaled_data = data_sections[1]
    scaling_func = data_sections[2]
    size(scaled_data, 2) >= 3 || throw(ArgumentError("scaled_data must have at least 3 columns: X, Y, E"))
    size(scaling_func, 2) >= 2 || throw(ArgumentError("scaling_func must have at least 2 columns: X_func, mu_func"))

    # Select the scaled scatter data columns
    X = Float64.(scaled_data[:, 1])
    Y = Float64.(scaled_data[:, 2])
    E = Float64.(scaled_data[:, 3]) # Y error
    all(isfinite, X) || throw(ArgumentError("scaled_data contains non-finite X values"))
    all(isfinite, Y) || throw(ArgumentError("scaled_data contains non-finite Y values"))
    all(e -> isfinite(e) && e > 0, E) || throw(ArgumentError("scaled_data contains non-finite or non-positive E values (cannot compute χ²)"))

    # Select the scaling function columns
    X_func = Float64.(scaling_func[:, 1]) # grid points of the scaling function
    mu_func = Float64.(scaling_func[:, 2]) # regression function value
    all(isfinite, X_func) || throw(ArgumentError("scaling_func contains non-finite X_func values"))
    all(isfinite, mu_func) || throw(ArgumentError("scaling_func contains non-finite mu_func values"))

    # Sort the scaling function data for later interpolation
    sort_idx = sortperm(X_func)
    X_sorted = X_func[sort_idx]
    mu_sorted = mu_func[sort_idx]

    # Filter out duplicate X points -- for interpolation and numerical derivatives
    X_unique = Float64[]
    mu_unique = Float64[]
    for i in eachindex(X_sorted)
        if isempty(X_unique) || X_sorted[i] != X_unique[end]
            push!(X_unique, X_sorted[i])
            push!(mu_unique, mu_sorted[i])
        end
    end
    length(X_unique) >= 2 || throw(ArgumentError("scaling_func must have at least 2 distinct X points for interpolation"))

    # Interpolate the scaling function at the data points
    itp = LinearInterpolation(X_unique, mu_unique, extrapolation_bc=Line())
    mu_at_data = itp.(X) # regression function value at the data points

    # Select the xerr column
    form = get(metadata, "form", 0)
    base_cols = form == 1 ? 8 : 7
    xerr = nothing
    if size(scaled_data, 2) == base_cols + 1
        xerr = Float64.(scaled_data[:, end])
        all(x -> isfinite(x) && x >= 0, xerr) || throw(ArgumentError("scaled_data contains non-finite or negative xerr values"))
    end

    xerr === nothing && sigma_xy !== nothing && throw(ArgumentError("sigma_xy is provided but xerr column is not present in scaled_data"))

    # case 1: no xerr column
    if xerr === nothing
        return sum(((Y .- mu_at_data) ./ E) .^ 2)
    end

    # Prepare the covariance vector
    sigma_xy_vec = if sigma_xy !== nothing
        vec = sigma_xy isa Real ? fill(Float64(sigma_xy), length(X)) : Float64.(sigma_xy)
        length(vec) == length(X) || throw(ArgumentError("sigma_xy length mismatch: expected $(length(X)), got $(length(vec))"))
        all(isfinite, vec) || throw(ArgumentError("sigma_xy contains non-finite values"))
        vec
    else
        nothing
    end

    # slopes[i] is the slope of the interval [X_unique[i], X_unique[i+1]]
    slopes = (mu_unique[2:end] .- mu_unique[1:end-1]) ./ (X_unique[2:end] .- X_unique[1:end-1])

    function dfdx_at(x::Float64)
        x <= X_unique[1] && return slopes[1]
        x >= X_unique[end] && return slopes[end]

        j = searchsortedlast(X_unique, x)
        j >= length(X_unique) && return slopes[end]

        X_unique[j] == x && 1 < j < length(X_unique) && return 0.5 * (slopes[j-1] + slopes[j])
        return slopes[j]
    end

    dfdx = map(dfdx_at, X)
    # case 2: xerr column is present
    sigma2_eff = E .^ 2 .+ (dfdx .^ 2) .* (xerr .^ 2)
    # case 3: sigma_xy column is present
    sigma_xy_vec !== nothing && (sigma2_eff .-= 2 .* dfdx .* sigma_xy_vec)
    all(s2 -> isfinite(s2) && s2 > 0, sigma2_eff) || throw(ArgumentError("Encountered non-finite or non-positive σ_eff^2 while computing χ²"))

    return sum((Y .- mu_at_data) .^ 2 ./ sigma2_eff)
end

"""
    chi2red_interp(metadata, data_sections; sigma_xy=nothing) -> Float64

Compute the reduced χ² using `chi2_interp`:

    χ²_red = χ² / (n_points - n_freeparams)

If `metadata["n_points"]` is missing, falls back to `size(data_sections[1], 1)`.
If `metadata["n_freeparams"]` is missing, defaults to 0.
"""
function chi2red_interp(metadata::Dict{String,Any},
                        data_sections::Vector{<:AbstractMatrix};
                        sigma_xy::Union{AbstractVector,Real,Nothing}=nothing)
    n_points = get(metadata, "n_points", size(data_sections[1], 1))
    n_freeparams = get(metadata, "n_freeparams", 0)
    dof = n_points - n_freeparams
    dof > 0 || throw(ArgumentError("Invalid degrees of freedom: n_points=$n_points, n_freeparams=$n_freeparams"))
    return chi2_interp(metadata, data_sections; sigma_xy=sigma_xy) / dof
end

"""
    extract_parameter_dict(metadata)

Convert BSA output p[:] array to unambiguous physical parameter names.

This is the "stdout" bridge - maps BSA's output convention back to parameter names
that match the input convention (Tc, c1, c2, c3).

# Parameter Mapping
- `scaling_form=0`: p[0]=Tc, p[1]=c1, p[2]=c2, p[3-5]=theta
- `scaling_form=1`: p[0]=Tc, p[1]=c1, p[2]=c3, p[3]=c2, p[4-8]=theta

# Returns
Dict{String, Tuple{Float64, Float64}} with keys: "Tc", "c1", "c2", "c3" (if form=1), "theta0-4"
Each value is a tuple (value, stderr).

# Data Flow
```
BSA output (.op file)
  ↓ parse_bsa_output()
metadata {"p0", "p1", "p2", ..., "form": 0}  ← ambiguous: p[2] depends on form
  ↓ extract_parameter_dict(metadata)
parameter_dict {"Tc", "c1", "c2", "c3", ...}  ← unambiguous parameter names
```

# Note
The reverse mapping (param_name → p[i]) is provided by `get_param_to_pidx_mapping()`.
"""
function extract_parameter_dict(metadata::Dict{String,Any})
    params = Dict{String,Tuple{Float64,Float64}}()  # (value, stderr)
    
    # Extract scaling_form from metadata
    scaling_form = get(metadata, "form", 0)
    
    # Common parameters (always p[0] and p[1])
    if haskey(metadata, "p0")
        params["Tc"] = metadata["p0"]
    end
    
    if haskey(metadata, "p1")
        params["c1"] = metadata["p1"]
    end
    
    # c2 and c3 mapping depends on scaling_form
    if scaling_form == 0
        # Standard form: p[2] = c2, p[3-5] = theta
        if haskey(metadata, "p2")
            params["c2"] = metadata["p2"]
        end
        if haskey(metadata, "p3")
            params["theta0"] = metadata["p3"]
        end
        if haskey(metadata, "p4")
            params["theta1"] = metadata["p4"]
        end
        if haskey(metadata, "p5")
            params["theta2"] = metadata["p5"]
        end
        
    elseif scaling_form == 1
        # With correction: p[2] = c3 (ω), p[3] = c2, p[4-8] = theta
        if haskey(metadata, "p2")
            params["c3"] = metadata["p2"]
        end
        if haskey(metadata, "p3")
            params["c2"] = metadata["p3"]
        end
        if haskey(metadata, "p4")
            params["theta0"] = metadata["p4"]
        end
        if haskey(metadata, "p5")
            params["theta1"] = metadata["p5"]
        end
        if haskey(metadata, "p6")
            params["theta2"] = metadata["p6"]
        end
        if haskey(metadata, "p7")
            params["theta3"] = metadata["p7"]
        end
        if haskey(metadata, "p8")
            params["theta4"] = metadata["p8"]
        end
    else
        @error "Invalid scaling_form: $scaling_form (must be 0 or 1)"
    end
    
    return params
end

"""
    extract_physical_params(parameter_dict; critical_param_name="Uc", eta_type=:none)

Convert unambiguous BSA parameters to user-friendly physical quantities.

This is the final step in the output processing chain - converts technical parameter
names (Tc, c1, c2) to physics-meaningful quantities with user-chosen naming/interpretation.

# Input
`parameter_dict::Dict{String, Tuple{Float64,Float64}}` from `extract_parameter_dict()`
Keys: "Tc", "c1", "c2", "c3" (if scaling_form=1), "theta0-4", ...

# Arguments
- `critical_param_name`: User-chosen name for critical point: "Uc", "Tc", "Jc", etc.
- `eta_type`: How to interpret c2: :eta_phi (η_φ = -1-c2), :eta_psi (η_ψ = -c2), or :none

# Returns
Dict{String, Tuple{Float64,Float64}} with user-friendly physical quantity names

# Data Flow
```
parameter_dict {"Tc" => (val, err), "c1" => (val, err), "c2" => (val, err)}
  ↓ extract_physical_params(parameter_dict, critical_param_name, eta_type)
physical_params {"Uc"/"Tc"/"Jc" => (val, err), "nu" => (val, err), "eta_psi" => (val, err)}
```

# Physical Transformations
- Tc → user-specified name (Uc/Tc/Jc): identity mapping, just rename
- c1 (1/ν) → ν: compute ν = 1/c1 with error propagation
- c2 → η_ψ or η_φ: η_ψ = -c2, η_φ = -1-c2 (user choice)
- c3 (if present) → ω: identity, just rename

# Example
```julia
params = extract_parameter_dict(metadata, 0)
# {"Tc" => (3.83, 0.005), "c1" => (0.94, 0.04), "c2" => (0.0, 0.0)}

phys = extract_physical_params(params, critical_param_name="Uc", eta_type=:none)
# {"Uc" => (3.83, 0.005), "nu" => (1.06, 0.045), "c2" => (0.0, 0.0)}
```
"""
function extract_physical_params(parameter_dict::Dict{String,Tuple{Float64,Float64}}; 
                                critical_param_name::String="Uc",
                                eta_type::Symbol=:none)
    phys = Dict{String,Tuple{Float64,Float64}}()
    
    # Critical point: Tc → Uc/Tc/Jc (user-specified name)
    if haskey(parameter_dict, "Tc")
        phys[critical_param_name] = parameter_dict["Tc"]
    end
    
    # Critical exponent: c1 → ν = 1/c1
    if haskey(parameter_dict, "c1")
        c1_val, c1_err = parameter_dict["c1"]
        nu_val = 1.0 / c1_val
        nu_err = c1_err / (c1_val^2)  # Error propagation: d(1/x)/dx = -1/x²
        phys["nu"] = (nu_val, nu_err)
    end
    
    # Anomalous dimension: c2 → η (user chooses interpretation)
    if haskey(parameter_dict, "c2")
        c2_val, c2_err = parameter_dict["c2"]
        
        if eta_type == :eta_phi
            phys["eta_phi"] = (-1.0 - c2_val, c2_err)
        elseif eta_type == :eta_psi
            phys["eta_psi"] = (-c2_val, c2_err)
        elseif eta_type == :none
            # Only output raw c2 when no interpretation is requested
            phys["c2"] = (c2_val, c2_err)
        end
    end
    
    # Correction exponent: c3 → ω (only present if scaling_form=1)
    if haskey(parameter_dict, "c3")
        phys["omega"] = parameter_dict["c3"]
    end
    
    return phys
end

"""
    print_summary(metadata; critical_param_name="Uc", eta_type=:none)

Pretty-print BSA fitting results following the three-layer architecture.

# Three-layer Output
1. **Metadata layer**: chi2, n_points, goodness of fit
2. **Params layer**: Unambiguous parameters (Tc, c1, c2, c3, ...)
3. **Phys layer**: User-friendly physical quantities (Uc/Tc/Jc, ν, η)

# Arguments
- `metadata`: Raw BSA output from `parse_bsa_output()` (contains "form" field for scaling_form)
- `critical_param_name`: Display name for critical point (default: "Uc")
- `eta_type`: Interpretation of c2 (:none, :eta_psi, :eta_phi)

# Example
```julia
metadata, _ = parse_bsa_output("result.op")
print_summary(metadata, critical_param_name="Uc", eta_type=:none)
```
"""
function print_summary(metadata::Dict; 
                      critical_param_name::String="Uc",
                      eta_type::Symbol=:none)
    println("\n" * "="^60)
    println("BSA Fitting Results Summary")
    println("="^60)
    
    # Extract parameter_dict and physical_params
    parameter_dict = extract_parameter_dict(metadata)
    physical_params = extract_physical_params(parameter_dict; 
                                             critical_param_name=critical_param_name,
                                             eta_type=eta_type)
    
    # Layer 2: Unambiguous parameters (Tc, c1, c2, ...)
    println("\n[Parameters]")
    for key in sort(collect(keys(parameter_dict)))
        val, err = parameter_dict[key]
        if !isnan(err) && err > 0
            rel_err = abs(err / val) * 100
            @printf("  %s = %.6e ± %.6e  (rel. err: %.2f%%)\n", key, val, err, rel_err)
        else
            @printf("  %s = %.6e\n", key, val)
        end
    end
    
    # Layer 3: Physical quantities (Uc, ν, η)
    if !isempty(physical_params)
        println("\n[Physical Quantities]")
        # Define display order for common physical quantities
        display_order = [critical_param_name, "nu", "eta_phi", "eta_psi", "c2", "omega"]
        for key in display_order
            if haskey(physical_params, key)
                val, err = physical_params[key]
                if !isnan(err) && err > 0
                    @printf("  %s = %.6f ± %.6f\n", key, val, err)
                else
                    @printf("  %s = %.6f\n", key, val)
                end
            end
        end
        # Print any remaining physical quantities not in display_order
        for key in sort(collect(keys(physical_params)))
            if !(key in display_order)
                val, err = physical_params[key]
                if !isnan(err) && err > 0
                    @printf("  %s = %.6f ± %.6f\n", key, val, err)
                else
                    @printf("  %s = %.6f\n", key, val)
                end
            end
        end
    end
    
    # Layer 1: Goodness of fit (from metadata)
    if haskey(metadata, "chi2") && haskey(metadata, "n_points")
        chi2 = metadata["chi2"]
        n_points = metadata["n_points"]
        chi2_red = chi2 / (n_points - 6)
        println("\n[Goodness of Fit]")
        @printf("  χ² = %.2f\n", chi2)
        @printf("  χ²_reduced = %.3f\n", chi2_red)
        if chi2_red < 0.5
            println("  → Excellent fit (possibly overfit)")
        elseif chi2_red < 1.5
            println("  → Good fit")
        elseif chi2_red < 3
            println("  → Acceptable fit")
        else
            println("  → Poor fit")
        end
    end
    
    println("="^60 * "\n")
end

end
