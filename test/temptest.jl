#using Jumol
include("../src/Jumol.jl")
println(pwd())
include("helper.jl")
molstruc, molcalc = create_toy_struct("n_aff_deform", 2.3, 0.5)
Jumol.initialize_potential!(molcalc,13)
Deformer = Jumol.Aff_deform(molstruc)
Jumol.calc_all_pair_forces!(molcalc, calc_stress_tensor = true)
quiet = true
Borelax = Jumol.Box_relax(molcalc, Deformer, max_num_iter = 100, energy_tol = 1.0e-5)
Borelax.stress_tol = 1.0e-8

#Jumol.relax_biaxial!(Borelax, 5.0e-2, 5.0e-2, quiet = quiet)
#Jumol.update_distances!(molstruc)
#@test sum([sum(molstruc.atom_list[k].distances) for k in 1:length(molstruc.atom_list)]) ≈ 5.606450178877528e8 atol = 1.0e2;
