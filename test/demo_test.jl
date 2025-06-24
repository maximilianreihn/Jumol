using Test 
using Suppressor
using Jumol

function demo_module_tests(exclude_time_intensiv)
    @testset "001_distance_calc_benchmarks folder tests" begin
        @suppress begin
            include("../demo_run/001_distance_calc_benchmarks/01_distance_4_atoms_benchmark.jl");
        end;
    end;

    @testset "002_potential_file_runs folder tests" begin
        @suppress begin
            include("../demo_run/002_potential_file_runs/01_LJ_6_12/01_LJ_6_12.jl")
            include("../demo_run/002_potential_file_runs/02_Tersoff_C/01_tersoff_carbon.jl")
            include("../demo_run/002_potential_file_runs/03_Yukawa_2D/01_yukawa_2D_silica.jl")
            include("../demo_run/002_potential_file_runs/04_Yukawa_extended/01_yukawa_2D_silica_baloon.jl")
        end;
    end;

    @testset "003_min_runs folder tests" begin 
        @suppress begin
            if !exclude_time_intensiv
                include("../demo_run/003_min_runs/01_minimize_4_atoms_lj_sd_vs_cg.jl");
                include("../demo_run/003_min_runs/02_minimize_4_atoms_sd_vs_cg.jl");
                include("../demo_run/003_min_runs/03_minimize_lattice_lj_sd_vs_cg.jl");
                include("../demo_run/003_min_runs/04_minimize_2D_silica_sample_sd_versus_cg.jl");
            end
        end;
    end;

    @testset "004_generate_lattices folder tests" begin 
        @suppress begin
            if !exclude_time_intensiv
                include("../demo_run/004_generate_lattices/01_LJ_glass_gen.jl");
                include("../demo_run/004_generate_lattices/02_MC_bond_switch_dual_harmonic.jl");
                include("../demo_run/004_generate_lattices/03_generate_mc_bond_switch.jl");
            end;
        end;
    end;

    @testset "005_affine_deformation folder tests" begin 
        @suppress begin
            if !exclude_time_intensiv
                include("../demo_run/005_affine_deformation/01_affine_deformation_4_atoms_benchmark.jl");
                include("../demo_run/005_affine_deformation/02_affine_deformation_4_atoms_video.jl");
                include("../demo_run/005_affine_deformation/03_affine_deformation_2D_silica_benchmark_example.jl");
                include("../demo_run/005_affine_deformation/04_affine_deformation_2D_silica_video.jl");

            end;
        end;
    end;

    @testset "006_AQS_run folder tests" begin
        @suppress begin
            if !exclude_time_intensiv
                include("../demo_run/006_AQS_run/01_box_relax/01_full_relax.jl");
                include("../demo_run/006_AQS_run/02_simple_shear/01_AQS_simple_shear_LJ_crystal.jl");
                include("../demo_run/006_AQS_run/02_simple_shear/02_AQS_simple_shear_LJ_glass.jl");
                include("../demo_run/006_AQS_run/02_simple_shear/03_AQS_simple_shear_network_glass.jl");
                include("../demo_run/006_AQS_run/02_simple_shear/04_AQS_simple_shear_network_glass_large.jl");
                include("../demo_run/006_AQS_run/02_simple_shear/05_AQS_linear_guesser_speed_up.jl");
                include("../demo_run/006_AQS_run/02_simple_shear/06_AQS_parallel_method.jl");
                include("../demo_run/006_AQS_run/03_tension/01_AQS_tension_LJ_glass.jl");
                include("../demo_run/006_AQS_run/03_tension/02_AQS_simple_tension_LJ_glass.jl");
            end;
        end;
    end;

    @testset "007_non_affine_disp_fields folder tests" begin 
        @suppress begin
            include("../demo_run/007_non_affine_disp_fields/01_calc_non_aff_disp_field.jl")
        end;
    end;

    @testset "008_int_run folder tests" begin 
        @suppress begin
            if !exclude_time_intensiv
                include("../demo_run/008_int_run/01_nve_run.jl");
                include("../demo_run/008_int_run/02_nve_lattice_lj.jl");
                include("../demo_run/008_int_run/03_nvt_4_atoms_benchmark.jl");
                include("../demo_run/008_int_run/04_nvt_LJ_glass_melt.jl");
                include("../demo_run/008_int_run/05_nvt_LJ_glass_quench.jl");
                include("../demo_run/008_int_run/06_nvt_run_network_glass.jl");
            end;
        end;
    end;

    @testset "009_image_silica folder tests" begin 
        @suppress begin
            if !exclude_time_intensiv
                include("../demo_run/009_image_silica/01_relax_image_sample.jl");
                include("../demo_run/009_image_silica/02_prepare_circular_test_samples.jl");
                # include("../demo_run/009_image_silica/03_shear_circular_sample.jl") #external files missing
            end;
        end;
    end;
end;
