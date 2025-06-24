## module defintion of JuMol
######################################
######################################
######################################

module Jumol

    import Core.Array
    using DelimitedFiles
    using LinearAlgebra
    using StaticArrays
    using Plots
    using PrettyTables
    using Random
    using Distributions
    using Statistics
    using Distributed
    using LinearSolve
    using BenchmarkTools
    using Profile        
    using Base.Threads
    using LaTeXStrings
    using Printf

    include("../src/sys/sys.jl")
    include("../src/structure/structure.jl")
    include("../src/structure/distance.jl")
    include("../src/calc/calc.jl")
    include("../src/calc/units.jl")
    include("../src/calc/integrator.jl")
    include("../src/calc/minimiser.jl")
    include("../src/structure/modifier.jl")
    include("../src/deform/deform.jl")
    include("../src/strategies/box_relax.jl")
    include("../src/calc/hessian.jl")
    include("../src/strategies/gen_lattices.jl")
    include("../src/strategies/aqs_procedure.jl")
    include("../src/strategies/aqs_parallel_procedure.jl")
    include("../src/strategies/aqs_parallel_procedure_threads.jl")
    include("../src/deform/disp_field.jl")
    include("../src/sys/AQS_acc_dev_basics.jl")

    function Base.:+(a::MVector{3,Float64}, b::MVector{3,Float64}) #inplace operation of MVectors in order to get much more efficient 
        @inbounds @simd for i in eachindex(a)
            a[i] += b[i]
        end
        return a
    end
    function Base.:-(a::MVector{3,Float64}, b::MVector{3,Float64}) #inplace operation of MVectors in order to get much more efficient 
        @inbounds @simd for i in eachindex(a)
            a[i] -= b[i]
        end
        return a
    end
end