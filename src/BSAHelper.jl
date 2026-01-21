module BSAHelper

include("bsa_core.jl")
include("bsa_bootstrap.jl")

using .BSACore: BSAConfig, BSAParameters, run_bsa_analysis, parse_bsa_output,
    extract_parameter_dict, extract_physical_params, print_summary, resolve_bsa_binary

using .BSABootstrap: BSAProblem, BootstrapConfig, BootstrapResult,
    bootstrap_bsa_analysis, save_bootstrap_summary, success_rate,
    prepare_bootstrap_plot_data, extract_and_format_physical_params

export BSAConfig, BSAParameters, BSAProblem, BootstrapConfig, BootstrapResult
export run_bsa_analysis, parse_bsa_output, extract_parameter_dict, extract_physical_params, print_summary
export bootstrap_bsa_analysis, save_bootstrap_summary, success_rate
export prepare_bootstrap_plot_data, extract_and_format_physical_params
export resolve_bsa_binary

# Optional plotting module (implemented via package extensions to avoid hard deps on PyPlot)
include("bsa_plotting.jl")

end
