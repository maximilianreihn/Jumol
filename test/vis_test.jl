using Test

function vis_module_tests(exclude_time_intensiv)
    @testset "network_vis.jl file implementation" begin
        num_x = 15
        num_y = 18
        bond_l = 3.05
        nw = Jumol.Network(bond_l, "harmonic_dual")
        allowed_rings = [4,5,6,7,8,9,10]
        alpha_target = [-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2]
        #alpha_target = [0.3 for k in 1:7]
        sig = 1
        target_stat = [0.0399705727098,0.275646469643,0.401653587686,0.212793419493,0.0578914101397,0.0101351126893,0.00132376850128]
        Jumol.initialize_hexagonal_network!(nw, num_x, num_y, sig, ring_size_list_allowed=allowed_rings, alpha_target_list=alpha_target, target_stat=target_stat)
        @test nw.Visnetwork.bond_lw == 2.0
        """
        Testing is done to see if there is an error creating the Visnetwork
        """

        molcalc = create_toy_struct("calc", 2.0, 0.4)
        molstruc = molcalc.structure
        Plotter = Jumol.Vis2d(molstruc)
        Plotter.bfig = 1000
        Plotter.hfig = 500

        Jumol.plot_covalent_bonds(Plotter)
        fig = Jumol.plot_box(Plotter)
        
        @test string(fig.attr) == vis_fig_attr
    end
end