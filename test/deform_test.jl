using Test

function deform_module_tests(exclude_time_intensiv= false)
    @testset "deform.jl file implementation" begin
        molstruc, molcalc = create_toy_struct("deform", 2.0, 0.4)
        Jumol.initialize_potential!(molcalc, 33)
        deformer = Jumol.Aff_deform(molstruc)
        Jumol.set_affine_deform_shear!(deformer,4.0,5.0,3.0)
        @test [molstruc.box.l1, molstruc.box.l2, molstruc.box.l3] == [5.0, 6.4031242374328485, 11.575836902790225]
        @test molstruc.box.e1 == Vector([1.0 ,0.0, 0.0])
        @test molstruc.box.e2 == Vector([0.6246950475544243, 0.7808688094430304, 0.0])
        @test molstruc.box.e3 == Vector([0.43193421279068006, 0.25916052767440806, 0.8638684255813601])

        molstruc, molcalc = create_toy_struct("deform", 2.0, 0.4)
        Jumol.initialize_potential!(molcalc, 33)
        deformer = Jumol.Aff_deform(molstruc)
        Jumol.set_affine_deform_vol!(deformer, 1.0, 2.0, 3.14)
        @test [molstruc.box.l1, molstruc.box.l2, molstruc.box.l3] == [6.0, 7.0, 13.14]

        molstruc, molcalc = create_toy_struct("deform", 2.0, 0.4)
        Jumol.initialize_potential!(molcalc, 33)
        deformer = Jumol.Aff_deform(molstruc)
        Jumol.set_affine_deform_true_shear!(deformer, 2.0, center_x = 3.14, center_y = 1.0)
        @test [molstruc.box.l1, molstruc.box.l2, molstruc.box.l3] == [7.0, 3.5714285714285716, 10.0]

        molstruc, molcalc = create_toy_struct("deform", 2.0, 0.4)
        Jumol.initialize_potential!(molcalc, 33)
        deformer = Jumol.Aff_deform(molstruc)
        Jumol.set_tension_to_atom_group_2d!(deformer, 1.0,-1.0,1.5, 1.2)
        @test [molstruc.box.l1, molstruc.box.l2, molstruc.box.l3] == [6.5, 6.2, 10.0]
    end
    @testset "disp_field.jl file implementation" begin
        molstruc, molcalc = create_toy_struct("n_aff_deform", 2.4, 0.4)
        Jumol.initialize_potential!(molcalc, 33)
        molstruc_non = deepcopy(molstruc)
        deform = Jumol.Aff_deform(molstruc)
        Jumol.set_affine_deform_shear!(deform, 1., 2., 3.)
        Jumol.update_distances!(molstruc)
        n_aff_disp = Jumol.Disp_field(molstruc, molstruc_non)
        Jumol.calc_disp_fields!(n_aff_disp)
        @test sum(abs.(n_aff_disp.nonaff_field)) ≈ 2.3448690905647496e-14 atol = 1.0e-8
        index = Jumol.locate_max_displacement(n_aff_disp, 99.)
        @test index == 68
    end
    @testset "strain_tensor_field.jl file implementation" begin
        #TODO
        """
        It seems that this part of the code is not used in any parts, urging to add unit tests if the code is acutally used
        """
    end
end