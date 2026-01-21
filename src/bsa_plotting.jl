module BSAPlotting
export build_title_from_phys, plot_bsa_data_collapse
export get_ylabel_text, latexstring_smart
export plot_data_points!, plot_window_indicator!

function get_default_markers end
function extract_base_variable end
function format_critical_param_latex end
function build_title_from_phys end
function latexstring_smart end
function get_ylabel_text end
function plot_scaling_function! end
function plot_data_points! end
function plot_residuals! end
function plot_window_indicator! end
function add_chi2_text! end
function plot_bsa_data_collapse end
end