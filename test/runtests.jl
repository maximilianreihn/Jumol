using Test
using Jumol
include("helper.jl")
include("calc_test.jl")
include("deform_test.jl")
include("strategies_test.jl")
include("structure_test.jl")
include("vis_test.jl")
include("demo_test.jl")

exclude_time_intensiv = true

@testset "unit" begin
    @testset "src/calc/ submodule" begin
        calc_module_tests(exclude_time_intensiv);
    end;

    @testset "src/deform/ submodule" begin
        deform_module_tests(exclude_time_intensiv)
    end;

    @testset "src/strategies/ submodule" begin
        strategies_module_tests(exclude_time_intensiv);
    end;

    @testset "src/structure/ submodule" begin
        structure_module_tests(exclude_time_intensiv);
    end;

    @testset "src/vis/ submodule" begin
        vis_module_tests(exclude_time_intensiv);
    end;

    @testset "demo_run/ module" begin
        demo_module_tests(exclude_time_intensiv)
    end;

end;
