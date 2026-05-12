using Printf
using Interpolations
using PyPlot
using LaTeXStrings

import BSAHelper: BSACore
import BSAHelper.BSAPlotting: add_chi2_text!, build_title_from_phys, extract_base_variable,
    format_critical_param_latex, get_default_markers, get_ylabel_text, latexstring_smart,
    plot_bsa_data_collapse, plot_data_points!, plot_residuals!, plot_scaling_function!,
    plot_window_indicator!

## -------------------------------------------------------------------------- ##
##                            Helper Functions                                ##
## -------------------------------------------------------------------------- ##

"""
    get_default_markers()

Return default marker cycle for different L values.
"""
function get_default_markers()
    return ["o", "s", "^", "v", "D", "<", ">", "p", "*", "h", "H", "d"]
end

"""
    extract_base_variable(critical_param_name::String)

Extract base variable name from critical parameter name.
Assumes format "Xc" → "X" (e.g., "Uc" → "U", "Tc" → "T", "Jc" → "J").

# Returns
- Base variable name (removes trailing 'c')
"""
function extract_base_variable(critical_param_name::String)
    if critical_param_name == "betac"
        return "\\beta"
    end
    if endswith(critical_param_name, "c")
        return critical_param_name[1:end-1]
    else
        return critical_param_name
    end
end

"""
    format_critical_param_latex(critical_param_name::String)

Convert critical parameter name to LaTeX format.
Assumes format "Xc" → "X_c" (e.g., "Uc" → "U_c", "Tc" → "T_c").

# Returns
- LaTeX-formatted parameter name
"""
function format_critical_param_latex(critical_param_name::String)
    if endswith(critical_param_name, "c") && length(critical_param_name) > 1
        base = critical_param_name[1:end-1]
        if base == "beta"
            base = "\\beta"
        end
        return "$(base)_c"
    else
        return critical_param_name
    end
end

"""
    build_title_from_phys(phys_fmt; critical_param_name="Uc", eta_type=:none)

Build title string from formatted physical quantities dictionary.

# Arguments
- `phys_fmt`: Dict{String, NamedTuple} from `extract_and_format_physical_params` or `format_physical_params`, containing
  pre-computed `value_str` and `error_str` with proper significant digits
- `critical_param_name`: Physical name for critical point (default: "Uc")
- `eta_type`: Interpretation of c2 (:none, :eta_phi, :eta_psi)

Returns a formatted title like: 
- "Data Collapse: Uc = 3.846 ± 0.012, ν = 1.06 ± 0.04"
"""
function build_title_from_phys(phys_fmt::Dict;
                               critical_param_name::String="Uc",
                               eta_type::Symbol=:none)
    title_parts = []
    
    # Helper function to add a quantity to title
    function add_quantity!(key::String, latex_name::String)
        if haskey(phys_fmt, key)
            entry = phys_fmt[key]
            val = entry.value
            err = entry.error
            # Skip display conditions
            # "none" is used when X-axis is correlation ratio (no critical point)
            if key == critical_param_name && critical_param_name == "none"
                return
            end
            if key == "nu" && (val == 0.0 || isinf(val) || isnan(val))
                return
            end
            if key == "c2" && val == 0.0
                return
            end
            # Use pre-formatted strings
            if !isnan(err) && err > 0
                push!(title_parts, "\$$(latex_name) = $(entry.value_str) \\pm $(entry.error_str)\$")
            else
                push!(title_parts, @sprintf("\$%s = %.4f\$", latex_name, val))
            end
        end
    end
    
    # Add critical parameter
    param_latex = format_critical_param_latex(critical_param_name)
    add_quantity!(critical_param_name, param_latex)
    
    # Add nu
    add_quantity!("nu", "\\nu")
    
    # Add eta based on eta_type
    if eta_type == :eta_phi
        add_quantity!("eta_phi", "\\eta_\\phi")
    elseif eta_type == :eta_psi
        add_quantity!("eta_psi", "\\eta_\\psi")
    elseif eta_type == :none
        add_quantity!("c2", "c_2")
    end
    
    # Add omega (c3) if present (scaling_form=1)
    add_quantity!("omega", "\\omega")
    
    return isempty(title_parts) ? "Data Collapse" : "Data Collapse: " * join(title_parts, ", ")
end

"""
    latexstring_smart(s::AbstractString)

Return `s` if `s` is a LaTeXString, otherwise convert `s` to LaTeXString.
"""
function latexstring_smart(s::AbstractString)
    return s isa LaTeXString ? s : L"%$(s)"
end

"""
    get_ylabel_text(observable_label, eta_type; c2_value=nothing)

Return Y-axis label based on eta_type and optional c2 value.

When `eta_type == :none`:
- If `c2_value == 0.0` (or not provided): no L scaling, returns `A`
- If `c2_value ≠ 0.0`: shows L scaling, returns `A / L^{c_2}`

Note: `observable_label` can be a plain String or LaTeXString. 
"""
function get_ylabel_text(observable_label::AbstractString, eta_type::Symbol; 
                        c2_value::Union{Float64, Nothing}=nothing)
    observable_label_latex = latexstring_smart(observable_label)

    # construct L suffix based on eta_type
    if eta_type == :none
        # Check if c2 is exactly zero (no scaling needed)
        l_suffix = c2_value == 0.0 ? "" : L" / L^{c_2}"
    elseif eta_type == :eta_phi
        l_suffix = L" / L^{-(1+\eta_\phi)}"
    elseif eta_type == :eta_psi
        l_suffix = L" / L^{-\eta_\psi}"
    else
        l_suffix = L" / L^{c_2}"
    end

    return observable_label_latex * l_suffix
end

## -------------------------------------------------------------------------- ##
##                            Plot parts of figure                            ##
## -------------------------------------------------------------------------- ##

"""
    plot_scaling_function!(ax, X_func, mu_func, sigma_func; 
                          label="Inferred \$F(X)\$", show_confidence=false)

Plot scaling function curve with optional confidence interval on given axis.

# Arguments
- `show_confidence`: If true, plot ±σ confidence band (default: false)
  Note: For Bootstrap analyses, σ from single fit doesn't reflect parameter uncertainty,
  so confidence band is disabled by default.
"""
function plot_scaling_function!(ax, X_func, mu_func, sigma_func; 
                                label::String="Inferred \$F(X)\$",
                                show_confidence::Bool=false)
    ax.plot(X_func, mu_func, "k-", linewidth=2.5, label=L"%$(label)", zorder=10)
    if show_confidence
        ax.fill_between(X_func, mu_func .- sigma_func, mu_func .+ sigma_func,
                        color="gray", alpha=0.3, label=L"$\pm\sigma$ confidence")
    end
end

"""
    plot_data_points!(ax, X, Y, E, L_values; mycolor=nothing, mymarker=nothing, kwargs...)

Plot data points grouped by L values with error bars.
"""
function plot_data_points!(ax, X::AbstractVector, Y::AbstractVector, E::AbstractVector, 
                          L::AbstractVector, L_values::AbstractVector;
                          mycolor=nothing, mymarker=nothing,
                          markersize=8, markeredgewidth=1.0, capsize=3, elinewidth=1, alpha=0.8,
                          xerr::Union{AbstractVector,Nothing}=nothing)
    # Use default markers if not provided
    default_markers = get_default_markers()
    
    for (idx, L_val) in enumerate(L_values)
        mask = L .== L_val
        if sum(mask) == 0
            continue
        end
        
        color = mycolor !== nothing ? mycolor(idx) : "C$(idx-1)"
        marker = mymarker !== nothing ? mymarker(idx) : default_markers[mod1(idx, length(default_markers))]
        
        if xerr !== nothing
            ax.errorbar(
                X[mask], Y[mask],
                xerr=xerr[mask], yerr=E[mask],
                fmt=marker,
                label=L"$L=%$(Int(L_val))$",
                color=color,
                markersize=markersize,
                markeredgewidth=markeredgewidth,
                mfc="none",
                capsize=capsize,
                elinewidth=elinewidth,
                alpha=alpha
            )
        else
            ax.errorbar(
                X[mask], Y[mask], yerr=E[mask],
                fmt=marker,
                label=L"$L=%$(Int(L_val))$",
                color=color,
                markersize=markersize,
                markeredgewidth=markeredgewidth,
                mfc="none",
                capsize=capsize,
                elinewidth=elinewidth,
                alpha=alpha
            )
        end
    end
end

"""
    plot_residuals!(ax, X, Y, E, L, L_values, X_func, mu_func;
                   mycolor=nothing, mymarker=nothing, xerr=nothing)

Plot residuals: (Y - F(X)) / E, where F is interpolated from scaling function.
"""
function plot_residuals!(ax, X::AbstractVector, Y::AbstractVector, E::AbstractVector, 
                        L::AbstractVector, L_values::AbstractVector,
                        X_func::AbstractVector, mu_func::AbstractVector;
                        mycolor=nothing, mymarker=nothing, xerr::Union{AbstractVector,Nothing}=nothing)
    # Create interpolation function
    itp = LinearInterpolation(X_func, mu_func, extrapolation_bc=Line())

    for (idx, L_val) in enumerate(L_values)
        mask = L .== L_val
        if sum(mask) == 0
            continue
        end
        
        X_data = X[mask]
        Y_data = Y[mask]
        E_data = E[mask]
        
        # Compute residuals
        Y_func_interp = itp.(X_data)
        residuals = (Y_data .- Y_func_interp) ./ E_data
        
        color = mycolor !== nothing ? mycolor(idx) : "C$(idx-1)"
        marker = mymarker !== nothing ? mymarker(idx) : "o"
        
        if xerr !== nothing
            xerr_data = xerr[mask]
            ax.errorbar(
                X_data, residuals,
                xerr=xerr_data, yerr=ones(length(residuals)),
                fmt=marker,
                color=color,
                markersize=6,
                markeredgewidth=0.8,
                mfc="none",
                capsize=2,
                elinewidth=0.8,
                alpha=0.8
            )
        else
            ax.errorbar(
                X_data, residuals, yerr=ones(length(residuals)),
                fmt=marker,
                color=color,
                markersize=6,
                markeredgewidth=0.8,
                mfc="none",
                capsize=2,
                elinewidth=0.8,
                alpha=0.8
            )
        end
    end
    
    # Add reference lines
    ax.axhline(0, color="k", linestyle="-", linewidth=1.5)
    ax.axhline(2, color="r", linestyle="--", linewidth=1, alpha=0.5)
    ax.axhline(-2, color="r", linestyle="--", linewidth=1, alpha=0.5)
end

"""
    plot_window_indicator!(ax_indicator, L_values;
                          all_sizes=[21,24,27,30,33,36,42,48,54,60,66,72],
                          mycolor=nothing, mymarker=nothing, show_xerr=true)

Create an innovative fitting window indicator on the given axis.
Highlights the sizes used in the fit with corresponding colors and markers.
Displays schematic error bars as legend-style indicators.

# Arguments
- `ax_indicator`: Axis for the window indicator
- `L_values`: Array of L values actually used in the fit
- `all_sizes`: All possible size values to show on the indicator
- `mycolor`: Color function/map (if nothing, uses default colors)
- `mymarker`: Marker function/map (if nothing, uses default markers)
- `show_xerr`: Whether to show horizontal error bars (default: true).
  Set to false for observables like Correlation Ratio that only have vertical errors.
"""
function plot_window_indicator!(ax_indicator, L_values::AbstractVector;
                                all_sizes::Vector{Int}=[21,24,27,30,33,36,42,48,54,60,66,72],
                                mycolor=nothing, mymarker=nothing, show_xerr::Bool=true)
    # Set up the indicator axis
    xmin = minimum(all_sizes) - 3
    xmax_data = maximum(all_sizes)
    extra_right = 6.0  # Extra space on the right for arrow + L label
    xmax = xmax_data + extra_right
    ax_indicator.set_xlim(xmin, xmax)
    ax_indicator.set_ylim(0, 1.2)  # Extended range to show markers fully
    ax_indicator.set_yticks([])
    ax_indicator.set_xticks(all_sizes)
    
    # Remove all tick marks (both major and minor)
    ax_indicator.tick_params(axis="x", which="both", length=0)
    
    # Configure axis: hide left/right/top spines, and use a custom arrow for the x-axis
    ax_indicator.spines["left"].set_visible(false)
    ax_indicator.spines["right"].set_visible(false)
    ax_indicator.spines["top"].set_visible(false)
    ax_indicator.spines["bottom"].set_visible(false)  # We draw the axis line manually

    # Draw a single arrow that serves as both the x-axis line and its arrow head.
    # The arrow runs from the left margin to slightly left of the right margin so
    # that there is room for the "L" label.
    axis_start = xmin
    axis_end   = xmax_data + 0.7 * extra_right
    ax_indicator.annotate("", xy=(axis_end, 0.0), xytext=(axis_start, 0.0),
                         arrowprops=Dict("arrowstyle"=>"->", "lw"=>1.5, "color"=>"black"),
                         annotation_clip=false)

    # Place the L label a bit to the right of the arrow head, inside the reserved
    # right-margin region so it is not squeezed at the figure boundary.
    label_x = axis_end + 0.1 * extra_right
    ax_indicator.text(label_x, 0.0, L"L", fontsize=14,
                      va="center", ha="left")
    
    # Default color and marker functions
    default_markers = get_default_markers()
    
    # Create tick labels with colors based on whether they're in the fitting window
    tick_labels = []
    tick_colors = []
    
    for size in all_sizes
        push!(tick_labels, string(size))
        if size in L_values
            # In fitting window: use black (high contrast)
            push!(tick_colors, "black")
        else
            # Outside fitting window: use light gray (low contrast)
            push!(tick_colors, "#BBBBBB")
        end
    end
    
    # Set tick labels
    ax_indicator.set_xticklabels(tick_labels, fontsize=10)
    
    # Color individual tick labels
    for (tick, color) in zip(ax_indicator.get_xticklabels(), tick_colors)
        tick.set_color(color)
    end
    
    # Add markers with schematic error bars above tick labels for sizes in fitting window
    for (idx, L_val) in enumerate(L_values)
        if L_val in all_sizes
            # Get color and marker for this L value
            color = mycolor !== nothing ? mycolor(idx) : "C$(idx-1)"
            marker = mymarker !== nothing ? mymarker(idx) : default_markers[mod1(idx, length(default_markers))]
            
            # Add schematic error bars with marker (legend-style, fixed size)
            # - Horizontal error bar (X direction): shown only when show_xerr=true
            # - Vertical error bar (Y direction): always shown
            xerr_val = show_xerr ? 1.0 : nothing
            ax_indicator.errorbar([L_val], [0.6], xerr=xerr_val, yerr=0.25,
                                 fmt=marker, color=color, ecolor=color,
                                 markersize=6, markeredgewidth=1, mfc="none",
                                 elinewidth=1, capsize=2, capthick=1, 
                                 alpha=0.8, zorder=3)
        end
    end

end

"""
    add_chi2_text!(ax, metadata::Dict)

Render the χ²_reduced text box. Prefers `metadata["chi2_eff"]` (σ_X-aware χ² that
`prepare_bootstrap_plot_data` injects at bootstrap time, picking either `chi2_interp`
or `chi2red_m2R` based on `eta_type`); falls back to `metadata["chi2"]` (BSA raw,
y-error only) when the corrected value is absent — e.g. when loading pre-1.1 JLD2
files that predate the chi2_eff injection.
"""
function add_chi2_text!(ax, metadata::Dict)
    haskey(metadata, "n_points") && haskey(metadata, "n_freeparams") || return
    chi2 = get(metadata, "chi2_eff", get(metadata, "chi2", nothing))
    chi2 === nothing && return

    dof = max(1, metadata["n_points"] - metadata["n_freeparams"])
    chi2_red = chi2 / dof
    ax.text(0.02, 0.95, @sprintf("\$\\chi^2_{\\mathrm{red}} = %.2f\$", chi2_red),
            transform=ax.transAxes, fontsize=14,
            verticalalignment="top",
            bbox=Dict("boxstyle"=>"round", "facecolor"=>"wheat", "alpha"=>0.8))
end

## -------------------------------------------------------------------------- ##
##                            Main Plotting Functions                         ##
## -------------------------------------------------------------------------- ##

"""
    plot_bsa_data_collapse(metadata, data_sections, phys_fmt, figs_dir; 
                          critical_param_name="Uc", eta_type=:none, kwargs...)

Render BSA data collapse plot with residuals using formatted physical quantities.

# Arguments
- `metadata`: BSA output metadata from `BSACore.parse_bsa_output()`
- `data_sections`: Data sections from BSA output
- `phys_fmt`: Dict{String, NamedTuple} from `format_physical_quantities`, containing
  pre-computed `value_str` and `error_str` with proper significant digits
- `figs_dir`: Output directory for figures
- `plot_mode`: Plot mode - `:full` (default) or `:simple`
  * `:full`: Shows data points, scaling function, residuals, title, and χ² info
  * `:simple`: Shows only data points (for publication-ready figures)
- `critical_param_name`: Physical name for critical point (default: "Uc")
- `eta_type`: Interpretation of c2 (:none, :eta_phi, :eta_psi)
- `observable_label`: Label for observable (default: "A"), accepts String or LaTeXString
- `xlabel_custom`: Custom x-axis label, accepts String or LaTeXString
- `show_confidence`: Whether to show ±σ confidence band (default: false)
  Set to true if using single BSA fit with MC errors. For Bootstrap analyses,
  keep false since σ doesn't reflect parameter uncertainty.
- Other styling parameters: `tick_params`, `font_legend`, `mycolor`, `mymarker`
"""
function plot_bsa_data_collapse(metadata::Dict, data_sections::Vector, phys_fmt::Dict,
                                figs_dir::String;
                                save_prefix::String="bsa_fss",
                                plot_mode::Symbol=:full,
                                tick_params::Dict=Dict(),
                                font_legend=nothing,
                                mycolor=nothing,
                                mymarker=nothing,
                                critical_param_name::String="Uc",
                                eta_type::Symbol=:none,
                                observable_label::AbstractString="A",
                                xlabel_custom::Union{AbstractString,Nothing}=nothing,
                                show_confidence::Bool=false)

    ## ------------------- Extract data from `data_sections` ------------------- ##
    
    if length(data_sections) < 2
        @error "Not enough data sections found in BSA output"
        return
    end

    scaled_data = data_sections[1]
    scaling_func = data_sections[2]
    has_correction = get(metadata, "form", 0) == 1

    ## ------------------------- Plot scatters and lines ------------------------ ##
    
    # Validate plot_mode
    if plot_mode ∉ [:full, :simple]
        @error "Unknown plot_mode: $plot_mode. Use :full or :simple"
    end
    
    # Extract data columns (with flexible handling of optional X errors column)
    # Expected columns from BSA:
    #   - form=0 (no correction): [X, Y, E, L, T, A, dA] = 7 cols
    #   - form=1 (with correction): [X, Y, E, X2, L, T, A, dA] = 8 cols
    # Our modification: if BSAProblem has x_err_col, an extra column is appended
    #   - form=0 + X errors: 8 cols
    #   - form=1 + X errors: 9 cols
    
    n_cols = size(scaled_data, 2)
    xerr = nothing
    
    if has_correction
        # form=1: expect 8 base columns, optionally 9 with X errors
        if n_cols == 9
            X, Y, E, X2, L, T, A, dA, xerr_data = eachcol(scaled_data)
            xerr = xerr_data
        elseif n_cols == 8
            X, Y, E, X2, L, T, A, dA = eachcol(scaled_data)
        else
            @error "Unexpected number of columns for form=1" n_cols expected="8 or 9"
            return
        end
    else
        # form=0: expect 7 base columns, optionally 8 with X errors
        if n_cols == 8
            X, Y, E, L, T, A, dA, xerr_data = eachcol(scaled_data)
            xerr = xerr_data
        elseif n_cols == 7
            X, Y, E, L, T, A, dA = eachcol(scaled_data)
        else
            @error "Unexpected number of columns for form=0" n_cols expected="7 or 8"
            return
        end
    end
    
    L_values = sort(unique(L))
    X_func, mu_func, sigma_func = eachcol(scaling_func)
    
    if plot_mode == :full
        # Full mode: data + scaling function + residuals
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9),
                                       gridspec_kw=Dict("height_ratios" => [3, 1],
                                                       "hspace" => 0.1))
        
        # Plot data points
        plot_data_points!(ax1, X, Y, E, L, L_values; mycolor=mycolor, mymarker=mymarker, xerr=xerr)
        
        # Plot scaling function
        plot_scaling_function!(ax1, X_func, mu_func, sigma_func; show_confidence=show_confidence)
        
        # Plot residuals
        plot_residuals!(ax2, X, Y, E, L, L_values, X_func, mu_func; 
                       mycolor=mycolor, mymarker=mymarker, xerr=xerr)
    
    else  # plot_mode == :simple
        # Simple mode: only data points with window indicator
        fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(6, 6.5),
                                       gridspec_kw=Dict("height_ratios" => [1, 12],
                                                       "hspace" => 0.15))
        
        # Plot fitting window indicator on ax0 (schematic error bars)
        plot_window_indicator!(ax0, L_values; mycolor=mycolor, mymarker=mymarker)
        
        # Plot data points on ax1
        plot_data_points!(ax1, X, Y, E, L, L_values; mycolor=mycolor, mymarker=mymarker, xerr=xerr)
    end

    ## ---------------------------- labels and title ---------------------------- ##
    
    # Get c2 value for ylabel determination (from phys_fmt if available)
    c2_value = haskey(phys_fmt, "c2") ? phys_fmt["c2"].value : nothing
    
    # Set labels
    ax1.set_ylabel(get_ylabel_text(observable_label, eta_type; c2_value=c2_value), fontsize=19)
    
    if plot_mode == :full
        # Full mode: set all labels, title, and chi2
        ax2.set_ylabel(L"Residuals / $\sigma$", fontsize=16)
        
        # X-axis label for bottom subplot
        if xlabel_custom !== nothing
            ax2.set_xlabel(latexstring_smart(xlabel_custom), fontsize=19)
        else
            base_var = extract_base_variable(critical_param_name)
            param_latex = format_critical_param_latex(critical_param_name)
            xlabel_text = latexstring("\$($base_var - $param_latex) L^{1/\\nu}\$")
            ax2.set_xlabel(xlabel_text, fontsize=19)
        end
        
        # Set title from formatted physical quantities
        ax1.set_title(build_title_from_phys(phys_fmt; critical_param_name=critical_param_name, 
                                            eta_type=eta_type), fontsize=19)
        
        # Add chi2 info
        add_chi2_text!(ax2, metadata)
    
    else  # plot_mode == :simple
        # Simple mode: set xlabel and add fit results text box in upper left
        if xlabel_custom !== nothing
            ax1.set_xlabel(latexstring_smart(xlabel_custom), fontsize=19)
        else
            base_var = extract_base_variable(critical_param_name)
            param_latex = format_critical_param_latex(critical_param_name)
            xlabel_text = latexstring("\$($base_var - $param_latex) L^{1/\\nu}\$")
            ax1.set_xlabel(xlabel_text, fontsize=19)
        end
        
        # Add fit results text box in upper left corner (replaces legend)
        fit_text = build_title_from_phys(phys_fmt; critical_param_name=critical_param_name, 
                                         eta_type=eta_type)
        # Remove "Data Collapse: " prefix for cleaner display
        fit_text = replace(fit_text, "Data Collapse: " => "")
        fit_text = replace(fit_text, ", " => "\n")
        ax1.text(0.05, 0.95, fit_text,
                transform=ax1.transAxes, fontsize=19,
                verticalalignment="top", horizontalalignment="left",
                bbox=Dict("boxstyle"=>"round", "facecolor"=>"white", "alpha"=>0.9, "edgecolor"=>"gray"))
    end
    
    ## ---------------------------------- style --------------------------------- ##

    if plot_mode == :full
        # Full mode: show legend
        ax1.legend(prop=font_legend, loc="best", framealpha=0.9)
    end
    # else: Simple mode - no legend (window indicator serves as legend)
    
    ax1.grid(linestyle="-", linewidth=0.5, alpha=0.4)
    ax1.minorticks_on()
    
    if plot_mode == :full
        # Full mode: apply style to both axes
        ax2.grid(linestyle="-", linewidth=0.5, alpha=0.4)
        ax2.minorticks_on()
        ax2.set_ylim(-5, 5)
        
        if !isempty(tick_params)
            ax1.tick_params(;tick_params...)
            ax2.tick_params(;tick_params...)
        end
    else
        # Simple mode: only style ax1
        if !isempty(tick_params)
            ax1.tick_params(;tick_params...)
        end
    end

    ## ------------------------------- save figure ------------------------------- ##

    plt.tight_layout()
    
    # Add mode suffix to filename
    mode_suffix = plot_mode == :simple ? "_simple" : ""
    fig_path_png = joinpath(figs_dir, "$(save_prefix)_collapse$(mode_suffix).png")
    
    fig.savefig(fig_path_png, dpi=300, bbox_inches="tight")
    display(fig)
    println("\n✓ Data collapse figure saved (mode: $plot_mode): $fig_path_png")
end

"""
    plot_bsa_data_collapse(metadata, data_sections, figs_dir; 
                          critical_param_name="Uc", eta_type=:none, kwargs...)

Convenience overload for single BSA fit (non-Bootstrap).

This method generates `phys_fmt` directly from `metadata` using simple
formatting (fixed significant digits for MC errors), then delegates to
the 4-argument version. Use this when you have only a single BSA fit
and don't need Bootstrap error estimation.

# Example
```julia
metadata, data_sections = BSACore.parse_bsa_output("output.op")
plot_bsa_data_collapse(metadata, data_sections, figs_dir;
    critical_param_name = "Uc",
    eta_type = :none
)
```
"""
function plot_bsa_data_collapse(metadata::Dict, data_sections::Vector,
                                figs_dir::String;
                                critical_param_name::String="Uc",
                                eta_type::Symbol=:none,
                                kwargs...)
    # Generate phys_fmt from metadata (single BSA fit with MC errors)
    parameter_dict = BSACore.extract_parameter_dict(metadata)
    phys_num = BSACore.extract_physical_params(parameter_dict;
                                               critical_param_name=critical_param_name,
                                               eta_type=eta_type)
    
    # Build phys_fmt with simple formatting (fixed 2 significant digits for errors)
    phys_fmt = Dict{String,NamedTuple}()
    for (key, (val, err)) in phys_num
        if err > 0
            # Simple formatting: match error precision to value
            val_str = @sprintf("%.4f", val)
            err_str = @sprintf("%.4f", err)
        else
            val_str = @sprintf("%.6f", val)
            err_str = @sprintf("%.6f", err)
        end
        phys_fmt[key] = (value = val, error = err, value_str = val_str, error_str = err_str)
    end
    
    # Delegate to the 4-argument version
    return plot_bsa_data_collapse(metadata, data_sections, phys_fmt, figs_dir;
                                  critical_param_name=critical_param_name,
                                  eta_type=eta_type,
                                  kwargs...)
end
