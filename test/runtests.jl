using BSAHelper
using Test

@testset "BSAHelper.jl" begin
    @testset "resolve_bsa_binary" begin
        tmp_path, io = mktemp()
        close(io)
        try
            @test resolve_bsa_binary(tmp_path) == tmp_path
        finally
            isfile(tmp_path) && rm(tmp_path; force=true)
        end
    end

    @testset "parse and map scaling_form=0" begin
        content = """
        # Number of data points = 4
        # Number of free parameters = 6
        # Scaling form: 0
        # chi^2 = 1.0
        # p[0] = 3.0 0.1
        # p[1] = 1.0 0.05
        # p[2] = 0.2 0.02
        # p[3] = 1.0 0.0
        # p[4] = 1.0 0.0
        # p[5] = 1.0 0.0

        1.0 2.0 0.1 4 3.9 2.0 0.1
        """
        path, io = mktemp()
        write(io, content)
        close(io)
        try
            metadata, data_sections = parse_bsa_output(path)
            @test metadata["form"] == 0
            @test length(data_sections) == 1

            params = extract_parameter_dict(metadata)
            @test isapprox(params["Tc"][1], 3.0)
            @test isapprox(params["c1"][1], 1.0)
            @test isapprox(params["c2"][1], 0.2)

            phys = extract_physical_params(params; critical_param_name="Uc", eta_type=:eta_phi)
            @test isapprox(phys["Uc"][1], 3.0)
            @test isapprox(phys["nu"][1], 1.0)
            @test isapprox(phys["eta_phi"][1], -1.2)
        finally
            isfile(path) && rm(path; force=true)
        end
    end

    @testset "parse and map scaling_form=1" begin
        content = """
        # Number of data points = 4
        # Number of free parameters = 7
        # Scaling form: 1
        # chi^2 = 2.0
        # p[0] = 2.5 0.1
        # p[1] = 0.8 0.05
        # p[2] = 0.3 0.01
        # p[3] = -0.1 0.02
        # p[4] = 1.0 0.0
        # p[5] = 1.0 0.0
        # p[6] = 1.0 0.0
        # p[7] = 1.0 0.0
        # p[8] = 1.0 0.0

        0.0 0.5 0.1 0.2 4 2.5 0.5 0.1
        """
        path, io = mktemp()
        write(io, content)
        close(io)
        try
            metadata, data_sections = parse_bsa_output(path)
            @test metadata["form"] == 1
            @test length(data_sections) == 1

            params = extract_parameter_dict(metadata)
            @test isapprox(params["Tc"][1], 2.5)
            @test isapprox(params["c1"][1], 0.8)
            @test isapprox(params["c3"][1], 0.3)
            @test isapprox(params["c2"][1], -0.1)

            phys = extract_physical_params(params; critical_param_name="Uc", eta_type=:eta_psi)
            @test isapprox(phys["Uc"][1], 2.5)
            @test isapprox(phys["nu"][1], 1.25)
            @test isapprox(phys["eta_psi"][1], 0.1)
            @test isapprox(phys["omega"][1], 0.3)
        finally
            isfile(path) && rm(path; force=true)
        end
    end

    @testset "plotting extension (optional)" begin
        pyplot_path = Base.find_package("PyPlot")
        latex_path = Base.find_package("LaTeXStrings")
        if pyplot_path === nothing || latex_path === nothing
            @info "PyPlot/LaTeXStrings not available; skipping plotting extension tests"
            @test true
        else
            @eval using PyPlot
            @eval using LaTeXStrings
            @test isdefined(BSAHelper, :BSAPlotting)
            @test length(methods(BSAHelper.BSAPlotting.plot_bsa_data_collapse)) > 0
        end
    end
end
