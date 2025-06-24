"""
This file includes the unit tests for the calc submodule of JuMol
All functions in the modul are (at the least implicitly) tested with this file
It is possible that due to different Julia versions or other dependencies the float64 results may differ slightly
"""
#TODO Debatable if the results should be checke with == or with ≈ 
using Test

#when testing you can exclude certain methods with this flag which take a long time, those tests are marked within an if statement
function calc_module_tests(exclude_time_intensiv = false)
    ##############################################
    @testset "calc.jl file implementation" begin
        """
        first all potentials and their different methods are tested
        In later stages differnt potential functions are possible however only one is tested since all methods within the potentials are tested in this 
        section already
        """
        #TODO every result which is an even 0.0 or 1.0 should be looked at 
        # println("calc.jl file implementation")
        potential_types = Dict("Harm_ring_ring"=>7,
            "Yukawa"=>13, 
            "Tersoff"=>21,
            "LJ_MF"=>33)

        for pot in keys(potential_types)
            if pot == "Harm_ring_ring"
                @testset "Harm ring ring potential" begin
                    calc_struct = create_toy_struct("calc")
                    force = zeros(Float64, 3)
                    @test Jumol.initialize_potential!(calc_struct, potential_types[pot]) ==  zeros(3,3)
                    @test calc_struct.potential.R_0_mat == hrr_ro_mat
                    @test Jumol.set_atom_forces_to_zero!(calc_struct) == zeros(3,3)
                    @test Jumol.calc_total_pot_en!(calc_struct) ≈ 14.222912360003367 atol = 1.0e-5
                    @test Jumol.calc_V(calc_struct.potential, 0.1,0.05) == 0.05
                    @test Jumol.calc_V_prime(calc_struct.potential, 0.1, 0.05) == 2.0
                    Jumol.calc_force!(force,calc_struct, 2, 4)
                    @test force == harmring_force
                    @test Jumol.calc_pot(calc_struct.potential, 5,5, 1.) ≈ 0.10292444847653455 atol = 1.0e-6
                    @test Jumol.calc_total_pot_en!(calc_struct) ≈ 14.222912360003367 atol = 1.0e-6
                    Jumol.calc_all_pair_forces!(calc_struct);
                    @test [calc_struct.structure.atom_list[k].force for k in 1:length(calc_struct.structure.atom_list)] == harmring_force_all
                    Jumol.calc_stress_tensor!(calc_struct);
                    @test calc_struct.stress_tensor == stress_harmring_a
                    
                end
            elseif pot == "Yukawa"
                @testset "Yukawa potential" begin
                    calc_struct = create_toy_struct("calc")
                    force = zeros(Float64, 3)
                    @test Jumol.initialize_potential!(calc_struct, potential_types[pot]) ==  zeros(3,3)
                    @test calc_struct.potential.V_prime_Rc == yukawa_vect
                    @test calc_struct.potential.params == yukawa_params
                    @test Jumol.calc_V_opt(calc_struct.potential, 1.2, 1.3, 3.14) ≈ 4.728922618137245 atol = 1.0e-6
                    @test Jumol.calc_V_prime_opt(calc_struct.potential, 1.2, 1.3, 3.14, 0.2) ≈ -3.338878405595603e8 atol = 1.0e-6
                    @test Jumol.calc_V_prime(calc_struct.potential, 1.2, 1.3, 3.14) ≈ -28.268151968287 atol = 1.0e-6
                    @test Jumol.calc_V(calc_struct.potential, 1.2, 1.3, 3.14) ≈ 4.728922618137247 atol = 1.0e-6
                    Jumol.calc_force!(force, calc_struct, 1, 2) 
                    @test force == yukawa_force
                    @test Jumol.calc_pot(calc_struct.potential, 1,2,calc_struct.potential.Rc/3) ≈ 57.18042090875043 atol = 1.0e-6
                    @test Jumol.calc_total_pot_en!(calc_struct) ≈ 4.732082467988773 atol=1.0e-6
                    Jumol.calc_all_pair_forces!(calc_struct);
                    @test [calc_struct.structure.atom_list[k].force for k in 1:length(calc_struct.structure.atom_list)] == forces_all
                    Jumol.calc_stress_tensor!(calc_struct);
                    @test calc_struct.stress_tensor == stress_tensor
                end
            elseif pot == "Tersoff"
                @testset "Tersoff potential" begin 
                    calc_struct = create_toy_struct("calc")
                    @test Jumol.initialize_potential!(calc_struct, potential_types[pot]) ==  zeros(3,3)
                    @test Jumol.calc_pot(calc_struct.potential, 1, 2) == -0.0
                    @test Jumol.calc_fc(0.1,calc_struct.potential.params) == 1
                    @test Jumol.calc_fc_fcp(0.1,calc_struct.potential.params) == (1.0, 0.0)
                    @test Jumol.calc_fr(0.1,calc_struct.potential.params) ≈ 983.2423255535989 atol = 1.0e-6
                    @test Jumol.calc_fr_frp(0.1,calc_struct.potential.params) == (983.2423255535989, -3429.4509072983974)
                    @test Jumol.calc_fa(0.1,calc_struct.potential.params) ≈ 277.93442916524674 atol = 1.0e-6
                    @test Jumol.calc_fa_fap(0.1,calc_struct.potential.params) == (277.93442916524674, -614.7631638706092)
                    @test Jumol.calc_b(calc_struct.potential,1,2) == 1.0
                    @test Jumol.calc_b_bp(calc_struct.potential,1,2) == (1.0, 0.0)
                    @test Jumol.calc_g_ijk(calc_struct.potential,1,2,3) == (0.0, 1.2959495026272982e6)
                    #TODO precomputed is not initalized @test Jumol.calc_dU_i_dr_ij(calc_struct.potential, 1, 2, 2.2, 1.2) == 1
                    #TODO precomputed is not initalized @test Jumol.calc_b_bp_full(calc_struct.potential, 1) == 1
                    @test Jumol.calc_total_pot_en!(calc_struct) ≈ -0.2372480456273478 atol=1.0e-6
                    Jumol.calc_all_pair_forces_Tersoff!(calc_struct);
                    @test [calc_struct.structure.atom_list[k].force for k in 1:length(calc_struct.structure.atom_list)] == tersoff_force_all
                end
            elseif pot == "LJ_MF"
                @testset "Lennard Jones potential" begin 
                    calc_struct = create_toy_struct("calc")
                    force = zeros(Float64, 3)
                    @test Jumol.initialize_potential!(calc_struct, potential_types[pot]) ==  zeros(3,3)
                    @test calc_struct.potential.V_prime_Rc == lj_prime_rc
                    @test calc_struct.potential.V_Rc == lj_rc
                    @test Jumol.calc_pot(calc_struct.potential, 1, 2, 0.001) ≈ 4.0e36 atol = 1.0e4
                    Jumol.calc_force!(force, calc_struct, 1,2) 
                    @test force == lj_force
                    @test Jumol.find_pos(calc_struct.potential, 1,2) == 2
                    @test Jumol.calc_total_pot_en!(calc_struct) ≈ -0.0009024729967356368 atol=1.0e-6
                    Jumol.calc_all_pair_forces!(calc_struct);
                    @test [calc_struct.structure.atom_list[k].force for k in 1:length(calc_struct.structure.atom_list)] == forces_all_b
                    Jumol.calc_stress_tensor!(calc_struct);
                    @test calc_struct.stress_tensor == stress_tensor_b
                end
            end
        end
    
    end
    ##############################################
    @testset "hessian.jl file implementation" begin
        # println("hessian.jl file implementation")
        molstruc, molcalc = create_toy_struct("hessian", 10.0, 0.5)
        Jumol.initialize_potential!(molcalc, 13)
        molhess = Jumol.Hessian(molcalc)
        ms2 = deepcopy(molstruc)
        deform = Jumol.Aff_deform(molstruc)
        Jumol.set_affine_deform_shear!(deform, 1., 1., 0.)
        Jumol.update_distances!(molstruc)
        Disp = Jumol.Disp_field(molstruc, ms2)
        Jumol.calc_disp_fields!(Disp)
        Jumol.calc_all_pair_forces!(molcalc, calc_stress_tensor = true)
        Jumol.calc_hessian_anal!(molhess)
        @test sum(abs.(molhess.hessian_mat)) == hessian_mat_sum
        if !exclude_time_intensiv
            Jumol.solve_eig!(molhess,calc_eigen_vector = true)
            @test sum(abs.(molhess.evals)) ≈ hessian_eig_sum atol = 1.0e-5
        end
    end
    ##############################################
    @testset "int.jl file implementation" begin
        # println("int.jl file implementation")
        molstruc, molcalc = create_toy_struct("integration", 2.4, 0.4)
        Jumol.initialize_potential!(molcalc, 33)
        molunits = Jumol.Units(unit_type=1)
        Jumol.set_units!(molunits)
        molint = Jumol.Integrator(molcalc,molunits, delta_t = 0.001)
        if !exclude_time_intensiv
            @test Jumol.run_nve!(molint, 1000, quiet = true) ≈ velocities  atol=1.0e-5
        end
    end
    ##############################################
    @testset "minimiser.jl file implementation" begin
        # println("minimiser.jl file implementation")
        alpha_min = 1.0e-1
        eps = 1.0e-7
        max_num_steps = 100000

        molstruc, molcalc = create_toy_struct("minimiser")
        Jumol.initialize_potential!(molcalc, 33)
        Minimizer = Jumol.Minimiser(molcalc)
        Minimizer.alpha_min = alpha_min
        Minimizer.acc_factor = 1.001
        Jumol.run_cg!(Minimizer,num_steps= max_num_steps, tolerance = eps, quiet = true)
        @test Minimizer.pot_en_hist[Minimizer.num_steps_real] ≈ -117.64866090673318 atol = 1.0e-6
        @test Minimizer.force_hist[Minimizer.num_steps_real] ≈ 8.597669399294241e-8 atol=1.0e-7

        molstruc, molcalc = create_toy_struct("minimiser")
        Jumol.initialize_potential!(molcalc, 33)
        Minimizer = Jumol.Minimiser(molcalc)
        Minimizer.alpha_min = alpha_min
        Jumol.run_sd!(Minimizer, max_num_steps, tolerance = eps, quiet = true)
        @test Minimizer.pot_en_hist[Minimizer.num_steps_real] ≈ 0.0 atol=1.0e-5
        @test Minimizer.force_hist[Minimizer.num_steps_real] ≈ 5.606186882196918e-13 atol = 1.0e-10

    end
end
