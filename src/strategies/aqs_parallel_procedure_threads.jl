
#########################################################
function aqs_sip_procedure_threads(shifter::Aqs_parallel_shifter)
    shifter.track_progress = 0
    starter_all = time()
    time_of_parall = 0.
    U_tracker = Jumol.U_obj[]
    starter_all_list = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, shifter.P)
    full_coordinate_history = Matrix{Float64}(undef, 3*shifter.minier.calc.structure.noa, shifter.aqs_steps)
    data_pool = Minimiser[deepcopy(shifter.minier) for _ in 1:shifter.P]
    sp_pool = Tuple{Vector{Float64}, Vector{Float64}}[(zeros(Float64,3*shifter.minier.calc.structure.noa),zeros(Float64, 9)) for k in 1:shifter.P]
    while shifter.track_progress < shifter.aqs_steps
        #println("Progress t=", shifter.track_progress, " time ", time()-starter_all, "s")
        coarse_range_left = Int64(ceil(minimum([shifter.P, (shifter.aqs_steps - shifter.track_progress)/shifter.D_split])))
        if coarse_range_left == 1
            last_steps = shifter.aqs_steps - shifter.track_progress
            start_last = time()
            inner_procedure!(q_num=shifter.q_num, num_steps=last_steps, 
                            eps= shifter.eps_mini, delta_gam=shifter.delta_gam, quiet=shifter.quiet, max_num_steps=shifter.max_min_steps, 
                            minier=shifter.minier, potential = shifter.potential, q_guesser_cache = shifter.q_cache)
            full_coordinate_history[:, shifter.track_progress+1:shifter.track_progress+last_steps] .= shifter.minier.rt_full[:,1:last_steps]
            
            time_of_parall += time()-start_last
            for final_acc in 1:last_steps
                @views u_temp = Jumol.U_obj(shifter.minier.U_hist[final_acc], norm(shifter.minier.rt_full[:, final_acc]), final_acc, 1, shifter.track_progress+final_acc,
                                    (1, shifter.track_progress), shifter.minier.step_hist[final_acc], shifter.minier.stress_hist[final_acc],  accepted = true)
                push!(U_tracker, u_temp)
            end
            shifter.track_progress += last_steps
        else
            starter_all_list[1] = (Jumol.get_global_atom_pos_vec(shifter.minier.calc.structure), Jumol.get_box_config(shifter.minier.calc.structure))
            for k in 2:coarse_range_left
                Jumol.set_affine_deform_shear!(shifter.minier.calc.structure, shifter.coarse_gam, 0.0, 0.0)
                Jumol.update_distances!(shifter.minier.calc.structure)
                starter_all_list[k] = (Jumol.get_global_atom_pos_vec(shifter.minier.calc.structure), Jumol.get_box_config(shifter.minier.calc.structure))
            end
            before_para = time()
            #parallel part here 
            #ddd = @allocated beg
                map(k->set_coordinates!(data_pool[k].calc.structure, (Jumol.get_global_atom_pos_vec(shifter.minier.calc.structure), Jumol.get_box_config(shifter.minier.calc.structure))), 1:shifter.P)
                @threads for coarse_step in 1:coarse_range_left
                    #iii = @allocated begin
                        st = time()
                        if coarse_step != 1
                            inner_cache = nothing
                        else
                            inner_cache = shifter.q_cache
                        end
                        Jumol.inner_procedure!(q_num=shifter.q_num, num_steps=shifter.D_split + Int64(coarse_step!=1), coarse_cacher = false,
                                    eps= shifter.eps_mini, delta_gam=shifter.delta_gam, quiet=shifter.quiet, max_num_steps=shifter.max_min_steps, 
                                    minier=data_pool[coarse_step], potential=shifter.potential, q_guesser_cache=inner_cache, is_inner = coarse_step!=1,
                                    init_position=starter_all_list[coarse_step])
                        sp_pool[coarse_step][1] .= Jumol.get_global_atom_pos_vec(data_pool[coarse_step].calc.structure)
                        sp_pool[coarse_step][2] .= Jumol.get_box_config(data_pool[coarse_step].calc.structure)
                    #end
                    #println("Alloc inside ", iii)
                end
            #end
            #println("The allocated in para was ", ddd)
            time_of_parall += time()-before_para
            #println("Para took ", time()-before_para)
            choosing_time = time()

            # choosing sequence starts here
            for check_index in 1:coarse_range_left-1
                if check_index==1
                    ending = shifter.D_split
                else
                    ending = shifter.D_split + 1
                end
                if abs(data_pool[check_index].U_hist[ending] - data_pool[check_index+1].U_hist[1]) > shifter.eps_accepter || shifter.track_progress + check_index*shifter.D_split >= shifter.aqs_steps
                    set_coordinates!(shifter.minier.calc.structure, sp_pool[check_index])
                    shifter.q_cache = nothing #because linear guessing is reset when plastic event happens
                    Jumol.push_dict_to_U_threads!(U_tracker, data_pool, shifter.track_progress, shifter.D_split, check_index, shifter.P, shifter.D_split)
                    #update_rts_hist!(full_coordinate_history, shifter, info_dict, check_index)
                    shifter.track_progress += check_index * shifter.D_split
                    break
                elseif check_index+1 == coarse_range_left
                    set_coordinates!(shifter.minier.calc.structure, sp_pool[check_index+1]) 
                    shifter.q_cache = @views Matrix{Float64}[data_pool[check_index+1].r_full[:, shifter.D_split-shifter.q_num+1:shifter.D_split], data_pool[check_index+1].rt_full[:, shifter.D_split-shifter.q_num+1:shifter.D_split]]
                    Jumol.push_dict_to_U_threads!(U_tracker, data_pool, shifter.track_progress, shifter.D_split, check_index+1, shifter.P, shifter.D_split)
                    #update_rts_hist!(full_coordinate_history, shifter, info_dict, check_index)
                    shifter.track_progress += (check_index+1) * shifter.D_split
                end
            end
        end
        if shifter.visualize
            plotter = Jumol.Vis2d(shifter.minier.calc.structure)
            fig = Jumol.plot_box(plotter)
            Jumol.plot_atomic_structure_binary(plotter, 10., 14.,
                                  [-shifter.minier.calc.structure.box.lx*0.075, shifter.minier.calc.structure.box.lx*1.075, shifter.minier.calc.structure.box.lx*1.075,-shifter.minier.calc.structure.box.lx*0.075],
                                  [-shifter.minier.calc.structure.box.ly*0.075, -shifter.minier.calc.structure.box.ly*0.075, shifter.minier.calc.structure.box.ly*1.075, shifter.minier.calc.structure.box.ly*1.075] )
            display(fig)
        end
    end
    
    sort!(U_tracker, by = x -> x.abs_index)   
    final_time = time()-starter_all
    overhead_time = final_time - time_of_parall
    return U_tracker, shifter.minier.calc.structure, full_coordinate_history, final_time, overhead_time
end

function push_dict_to_U_threads!(U_tracker::Vector{Jumol.U_obj}, data_pool::Vector{Jumol.Minimiser}, absolut_tracker::Int64, fine_steps::Int64, boarder, PP::Int64, D_split::Int64)
    """
    gets data of all parallel aqs steps from the info dict and creates a u_obj for every aqs step and marks them as accepted or denied
    """;
    for coarse_k in 1:PP
        if coarse_k<=boarder
            status = true
        else
            status = false
        end
        if coarse_k==1
            ending = D_split
        else
            ending = D_split + 1
        end
        giver = (coarse_k-1, absolut_tracker) #-1 in order to make it comform with plots
        for k in 1:ending
            if coarse_k > 1
                temp_absolut_tracker = absolut_tracker - 1
            else
                temp_absolut_tracker = absolut_tracker
            end
            #println("This is ",k, " ", coarse_k, " ", giver, " ", status)
            @views begin
                push!(U_tracker, U_obj(data_pool[coarse_k].U_hist[k],
                                norm(data_pool[coarse_k].rt_full[:,k]),
                                k, 
                                coarse_k,
                                temp_absolut_tracker + fine_steps*(coarse_k-1)+ k,
                                giver,
                                data_pool[coarse_k].step_hist[k],
                                data_pool[coarse_k].stress_hist[k],
                                accepted = status))
            end
        end
    end
end
