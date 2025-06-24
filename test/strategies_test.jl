
using Test
using Random

Random.seed!(999)

function structure_module_tests(exclude_time_intensiv)
    @testset "network.jl file implementation" begin
        num_x = 25
        num_y = 26
        bond_l = 3.05
        nw = Jumol.Network(bond_l, "harmonic_dual")
        allowed_rings = [4,5,6,7,8,9,10]
        alpha_target = [-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2]
        #alpha_target = [0.3 for k in 1:7]
        sig = 1
        target_stat = [0.0399705727098,0.275646469643,0.401653587686,0.212793419493,0.0578914101397,0.0101351126893,0.00132376850128]
        Jumol.initialize_hexagonal_network!(nw, num_x, num_y, sig, ring_size_list_allowed=allowed_rings, alpha_target_list=alpha_target, target_stat=target_stat)
        #if !exclude_time_intensiv
        if true
            Jumol.run_dual_MC_bond_switch(25, nw, num_x, num_y, plot_it = false, save = false)
            @test sum(nw.Ring_stat.ring_stat_mat) == 29.0
        end
        
        """
        Jumol.switch_bond_set seems to not be used in any code, if usage appears, add unit test as well
        """
    end
    @testset "Gen_lattice.jl file implementation" begin
        molstruc = Jumol.Structure()
        molstruc.rc = 2.4;molstruc.rskin = 0.4;molstruc.pbx = molstruc.pby= 1
        
        Generate_lat = Jumol.Gen_lattice(molstruc)
        Jumol.create_LJ_lattice!(Generate_lat, 6, 3)
        @test [[molstruc.atom_list[k].pos, molstruc.atom_list[k].vel,molstruc.atom_list[k].acc] for k in 1:length(molstruc.atom_list)] == lj_glass_atoms
        
        Random.seed!(999)
        molstruc = Jumol.Structure()
        molstruc.rc = 2.4;molstruc.rskin = 0.4;molstruc.pbx = molstruc.pby= 1
        
        Generate_lat = Jumol.Gen_lattice(molstruc)
        Jumol.create_random_pos_LJ!(Generate_lat, 6, 3)
        @test [[molstruc.atom_list[k].pos, molstruc.atom_list[k].vel,molstruc.atom_list[k].acc] for k in 1:length(molstruc.atom_list)] == lj_glass_atoms_random
    end
    @testset "box_relax.jl file implementation" begin
        molstruc, molcalc = create_toy_struct("n_aff_deform", 2.3, 0.5)
        Jumol.initialize_potential!(molcalc,13)
        Deformer = Jumol.Aff_deform(molstruc)
        Jumol.calc_all_pair_forces!(molcalc, calc_stress_tensor = true)
        @test sum([sum(molstruc.atom_list[k].distances) for k in 1:length(molstruc.atom_list)]) ≈ -4845.015583232061 atol = 1.0e-2
        quiet = true
        Borelax = Jumol.Box_relax(molcalc, Deformer, max_num_iter = 100, energy_tol = 1.0e-5)
        Borelax.stress_tol = 1.0e-8
        #if !exclude_time_intensiv
        if false # volume deformation broken, todo
            Jumol.relax_biaxial!(Borelax, 5.0e-2, 5.0e-2, quiet = quiet)
            Jumol.update_distances!(molstruc)
            @test sum([sum(molstruc.atom_list[k].distances) for k in 1:length(molstruc.atom_list)]) ≈ 5.606450178877528e8 atol = 1.0e2;
        end
    end
    @testset "aqs_procedure.jl file implementation" begin
        molstruc, molcalc = create_toy_struct("aqs_tester", 2.4, 0.4)
        Jumol.initialize_potential!(molcalc, 33)
        shifter_aqs = Jumol.Aqs_shifter(molcalc)
        results = Jumol.aqs_procedure!(shifter_aqs)
        @test results[2] ≈ -393.4861447979207  atol = 1.0e-4;
    end
end