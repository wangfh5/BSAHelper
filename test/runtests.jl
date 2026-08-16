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

    @testset "named point predictions (scaling_form=1)" begin
        columns = "row X1 X2 Y E mu_full std_full_latent std_full_legacy mu_zero std_zero_latent std_zero_legacy cov_full_zero_latent correction std_correction_latent dmu_full_dX1 dmu_zero_dX1 L T A dA"
        content = """
        # Scaling form : 1
        # Number of data points in dataset[0] = 2
        # Number of free parameters = 1

        0.0 1.0 0.1 0.2 4 0.0 1.0 0.1
        1.0 2.0 0.2 0.1 8 1.0 2.0 0.2


        0.0 0.9 0.1
        1.0 1.8 0.1


        0.0 1.0 0.1 0.2 4 0.0 1.0 0.1
        1.0 2.0 0.2 0.1 8 1.0 2.0 0.2


        0.0 0.9 0.1
        1.0 1.8 0.1


        # Dataset : 0
        # Section : point_predictions
        # Columns : $columns
        0 0.0 0.2 1.0 0.1 0.95 0.02 0.11 0.90 0.01 0.10 0.0001 0.05 0.02 1.0 0.8 4 0.0 1.0 0.1
        1 1.0 0.1 2.0 0.2 1.90 0.03 0.21 1.80 0.02 0.20 0.0002 0.10 0.03 0.5 0.4 8 1.0 2.0 0.2
        """
        path, io = mktemp()
        write(io, content)
        close(io)
        try
            metadata, data_sections = parse_bsa_output(path)
            @test metadata["section_index"][(0, "point_predictions")] == 5
            @test has_point_predictions(metadata)
            predictions = get_point_predictions(metadata, data_sections)
            @test predictions.row == [0.0, 1.0]
            @test predictions.correction ≈ predictions.mu_full - predictions.mu_zero
            @test chi2_interp(metadata, data_sections) ≈
                  sum(((predictions.Y - predictions.mu_full) ./ predictions.E) .^ 2)

            data_sections[1] = hcat(data_sections[1], [0.03, 0.04])
            sigma_xy = [0.0002, -0.0001]
            expected_variance = predictions.E .^ 2 .+
                (predictions.dmu_full_dX1 .* data_sections[1][:, end]) .^ 2 .-
                2 .* predictions.dmu_full_dX1 .* sigma_xy
            expected = sum((predictions.Y - predictions.mu_full) .^ 2 ./ expected_variance)
            @test chi2_interp(metadata, data_sections; sigma_xy=sigma_xy) ≈ expected
            @test chi2red_interp(metadata, data_sections; sigma_xy=sigma_xy) ≈ expected
        finally
            isfile(path) && rm(path; force=true)
        end
    end

    @testset "named form-0 point predictions" begin
        columns = "row X1 Y E mu std_latent dmu_dX1 L T A dA"
        header = """
        # Scaling form : 0
        # Number of data points in dataset[0] = 2
        # Number of free parameters = 1

        0.0 1.0 0.1 4 0.0 1.0 0.1
        1.0 2.0 0.2 8 1.0 2.0 0.2


        0.0 0.9 0.1
        1.0 1.8 0.1
        """
        predictions_block = """

        # Dataset : 0
        # Section : point_predictions
        # Columns : $columns
        0 0.0 1.0 0.1 0.95 0.02 1.0 4 0.0 1.0 0.1
        1 1.0 2.0 0.2 1.90 0.03 0.5 8 1.0 2.0 0.2
        """
        path, io = mktemp()
        write(io, header * predictions_block)
        close(io)
        try
            metadata, data_sections = parse_bsa_output(path)
            @test metadata["form"] == 0
            @test has_point_predictions(metadata)
            predictions = get_point_predictions(metadata, data_sections)
            @test predictions.mu == [0.95, 1.90]
            # χ² uses the native per-point predictions (mu), not the
            # interpolated grid values ([0.9, 1.8])
            @test chi2_interp(metadata, data_sections) ≈
                  sum(((predictions.Y - predictions.mu) ./ predictions.E) .^ 2)
            @test chi2_interp(metadata, data_sections) ≉ 2.0

            data_sections[1] = hcat(data_sections[1], [0.03, 0.04])
            sigma_xy = [0.0002, -0.0001]
            expected_variance = predictions.E .^ 2 .+
                (predictions.dmu_dX1 .* data_sections[1][:, end]) .^ 2 .-
                2 .* predictions.dmu_dX1 .* sigma_xy
            expected = sum((predictions.Y - predictions.mu) .^ 2 ./ expected_variance)
            @test chi2_interp(metadata, data_sections; sigma_xy=sigma_xy) ≈ expected
            @test chi2red_interp(metadata, data_sections; sigma_xy=sigma_xy) ≈ expected
        finally
            isfile(path) && rm(path; force=true)
        end

        # Legacy form-0 output without a point_predictions section falls back
        # to interpolating the scaling-function grid
        path, io = mktemp()
        write(io, header)
        close(io)
        try
            metadata, data_sections = parse_bsa_output(path)
            @test !has_point_predictions(metadata)
            @test chi2_interp(metadata, data_sections) ≈ 2.0
        finally
            isfile(path) && rm(path; force=true)
        end
    end

    @testset "legacy point_predictions_form1 alias" begin
        content = """
        # Scaling form : 1
        # Dataset : 0
        # Section : point_predictions_form1
        # Columns : row X1
        0 0.2
        1 0.3
        """
        path, io = mktemp()
        write(io, content)
        close(io)
        try
            metadata, data_sections = parse_bsa_output(path)
            @test metadata["section_index"][(0, "point_predictions_form1")] == 1
            @test has_point_predictions(metadata)
            @test get_point_predictions(metadata, data_sections).X1 == [0.2, 0.3]
        finally
            isfile(path) && rm(path; force=true)
        end
    end

    @testset "physical parameter formatting quotes order-of-magnitude error with warning" begin
        # err_of_err comparable to err (too few successful bootstrap samples):
        # the error is escalated to its order of magnitude and a warning is emitted
        formatted = @test_logs (:warn, r"order of magnitude") BSAHelper.BSABootstrap.format_physical_params(
            Dict("eta_phi" => (0.6691286447059479, 0.009151253245169017)), 10)
        @test formatted["eta_phi"].digits == 0
        @test formatted["eta_phi"].value_str == "0.67"
        @test formatted["eta_phi"].error_str == "0.01"

        # Integer-valued Uc: value precision stays aligned with the escalated
        # error instead of collapsing to "5" ± "0"
        formatted = @test_logs (:warn, r"order of magnitude") BSAHelper.BSABootstrap.format_physical_params(
            Dict("Uc" => (5.423756758475438, 0.0069434274377905714)), 10)
        @test formatted["Uc"].digits == 0
        @test formatted["Uc"].value_str == "5.42"
        @test formatted["Uc"].error_str == "0.01"
    end

    @testset "physical parameter formatting keeps significant error digits when well sampled" begin
        formatted = BSAHelper.BSABootstrap.format_physical_params(
            Dict("eta_phi" => (0.6691286447059479, 0.009151253245169017)), 1000)
        @test formatted["eta_phi"].digits == 1
        @test formatted["eta_phi"].value_str == "0.67"
        @test formatted["eta_phi"].error_str == "0.01"
    end
end
