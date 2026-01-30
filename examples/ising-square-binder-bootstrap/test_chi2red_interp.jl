using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using BSAHelper

## -------------------------------------------------------------------------- ##
##                              Setup and Data                                ##
## -------------------------------------------------------------------------- ##

const OP_PATH = joinpath(@__DIR__, "results", "bootstrap_recon.op")
isfile(OP_PATH) || error("Missing example output file: $OP_PATH")

metadata, data_sections = parse_bsa_output(OP_PATH)
scaled_data = data_sections[1]
n_rows = size(scaled_data, 1)

n_points = get(metadata, "n_points", n_rows)
n_freeparams = get(metadata, "n_freeparams", 0)
dof = n_points - n_freeparams

println("="^60)
println("Test: chi2_interp / chi2red_interp")
println("="^60)
println("Data file: $OP_PATH")
println("n_points = $n_points, n_freeparams = $n_freeparams, dof = $dof")
println()

## -------------------------------------------------------------------------- ##
##                        1. Basic chi2 calculation                           ##
## -------------------------------------------------------------------------- ##

println("[1] Basic chi2 calculation (no xerr)")

chi2 = chi2_interp(metadata, data_sections)
chi2red = chi2red_interp(metadata, data_sections)

@assert isfinite(chi2) && chi2 > 0 "chi2 must be finite and positive"
@assert isfinite(chi2red) && chi2red > 0 "chi2red must be finite and positive"
@assert dof > 0 "degrees of freedom must be positive"
@assert chi2red == chi2 / dof "chi2red must equal chi2 / dof"

# Regression check
const EXPECTED_CHI2 = 1890.9910095085772
@assert isapprox(chi2, EXPECTED_CHI2; atol=1e-6) "Regression test failed"

println("    chi2     = $chi2")
println("    chi2_red = $chi2red")
println("    OK")
println()

## -------------------------------------------------------------------------- ##
##                      2. X-error propagation (xerr)                         ##
## -------------------------------------------------------------------------- ##

println("[2] X-error propagation: σ_eff² = σ_y² + (df/dx)² σ_x²")

# xerr = 0 should give the same result as no xerr
ds_xerr0 = copy(data_sections)
ds_xerr0[1] = hcat(scaled_data, zeros(n_rows))
chi2_xerr0 = chi2_interp(metadata, ds_xerr0)
@assert isapprox(chi2_xerr0, chi2; atol=1e-12) "xerr=0 should match original chi2"
println("    xerr = 0: chi2 = $chi2_xerr0 (matches original)")

# xerr > 0 should reduce chi2 (larger effective variance)
ds_xerr = copy(data_sections)
ds_xerr[1] = hcat(scaled_data, fill(0.05, n_rows))
chi2_xerr = chi2_interp(metadata, ds_xerr)
chi2red_xerr = chi2red_interp(metadata, ds_xerr)

@assert isfinite(chi2_xerr) && chi2_xerr > 0
@assert chi2red_xerr == chi2_xerr / dof
@assert chi2_xerr < chi2 "xerr > 0 should reduce chi2"

println("    xerr = 0.05: chi2 = $chi2_xerr (reduced from $chi2)")
println("    chi2_red = $chi2red_xerr")
println("    OK")
println()

## -------------------------------------------------------------------------- ##
##                    3. Covariance term (sigma_xy)                           ##
## -------------------------------------------------------------------------- ##

println("[3] Covariance term: σ_eff² = σ_y² + (df/dx)² σ_x² - 2(df/dx) σ_xy")

# sigma_xy = 0 should give the same result as no sigma_xy
sigma_xy_zero = zeros(n_rows)
chi2_cov0 = chi2_interp(metadata, ds_xerr; sigma_xy=sigma_xy_zero)
@assert isapprox(chi2_cov0, chi2_xerr; atol=1e-12) "sigma_xy=0 should match xerr-only chi2"
println("    sigma_xy = 0 (vector): chi2 = $chi2_cov0 (matches xerr-only)")

# scalar sigma_xy = 0 should also work
chi2_cov_scalar0 = chi2_interp(metadata, ds_xerr; sigma_xy=0.0)
@assert isapprox(chi2_cov_scalar0, chi2_xerr; atol=1e-12) "scalar sigma_xy=0 should match xerr-only chi2"
println("    sigma_xy = 0 (scalar): chi2 = $chi2_cov_scalar0 (matches xerr-only)")
println("    OK")
println()

## -------------------------------------------------------------------------- ##
##                         4. Error handling                                  ##
## -------------------------------------------------------------------------- ##

println("[4] Error handling")

# sigma_xy length mismatch should throw ArgumentError
bad_sigma_xy = zeros(n_rows - 1)
try
    chi2_interp(metadata, ds_xerr; sigma_xy=bad_sigma_xy)
    @assert false "Should have thrown ArgumentError"
catch e
    @assert e isa ArgumentError "Expected ArgumentError, got $(typeof(e))"
    println("    sigma_xy length mismatch: ArgumentError thrown correctly")
end

# sigma_xy without xerr should throw ArgumentError
try
    chi2_interp(metadata, data_sections; sigma_xy=sigma_xy_zero)
    @assert false "Should have thrown ArgumentError"
catch e
    @assert e isa ArgumentError "Expected ArgumentError, got $(typeof(e))"
    println("    sigma_xy without xerr: ArgumentError thrown correctly")
end

println("    OK")
println()

## -------------------------------------------------------------------------- ##
##                              Summary                                       ##
## -------------------------------------------------------------------------- ##

println("="^60)
println("All tests passed!")
println("="^60)
