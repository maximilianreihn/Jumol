mutable struct Aqs_parallel_shifter
    minier::Minimiser
    aqs_steps::Int64
    max_min_steps::Int64
    eps_mini::Float64
    delta_gam::Float64
    dist_break::Union{Float64, Nothing}
    linear_guessing::Bool
    q_num::Int64
    quiet::Bool
    visualize::Bool
    siO_glass::Bool
    potential::Int64
    P::Int64
    D_split::Int64
    eps_accepter::Float64
    q_cache::Union{Nothing,Vector{Matrix{Float64}}}
    track_progress::Int64
    coarse_gam::Float64
    coarse_aqs_steps::Int64
    function Aqs_parallel_shifter(minier;
        aqs_steps = 30,
        max_min_steps=100000,
        eps_mini = 1.0e-6,
        delta_gam = 0.005,
        dist_break = 0.5,
        linear_guessing= true,
        q_num=0,
        quiet = true,
        visualize = false,
        siO_glass = false,
        potential = 33,
        P = 4,
        D_split = 10,
        eps_accepter = 2.0e-6,
        q_cache=nothing,
        track_progress=0)

        coarse_gam = D_split * delta_gam 
        coarse_aqs_steps = ceil(aqs_steps/D_split)

        return new(minier,
        aqs_steps, max_min_steps, eps_mini, delta_gam, dist_break,
        linear_guessing, q_num, quiet, visualize, siO_glass, potential, 
        P, D_split, eps_accepter, q_cache,track_progress, coarse_gam, coarse_aqs_steps)
    end
end

function get_box_config(structure::Structure)
    return Float64[structure.box.lx,structure.box.ly,structure.box.lz,
    structure.box.lxy,structure.box.lyz,structure.box.lxz,
    structure.box.l1,structure.box.l2,structure.box.l3]
end

function set_box_config!(structure::Structure, box_config::Vector{Float64})
    structure.box.lx,structure.box.ly,structure.box.lz,
    structure.box.lxy,structure.box.lyz,structure.box.lxz,
    structure.box.l1,structure.box.l2,structure.box.l3 = box_config
end

function best_approx_full(r_0::Vector{Float64}, r_i::Vector{Float64}, base_k::Matrix{Float64})
    return transpose(base_k) * (r_i-r_0)
end


function aqs_sequentiell!(;
    q_num::Int64=0,
    num_steps::Int64= 35,
    delta_gam::Float64 = 0.005,
    quiet::Bool= true,
    coarse_cacher::Bool= false,
    dist_break::Union{Nothing,Float64}= nothing,
    eps::Float64= 1.0e-5,
    max_num_steps::Int64= 100000,
    minier::Union{Nothing,Minimiser}= nothing,
    potential::Union{Nothing,Int64}= nothing,
    siO_glass::Bool= false,
    q_guesser_cache::Union{Nothing,Vector{Matrix{Float64}}}= nothing,
    mini_breaker_search::Union{Nothing,Float64}= nothing, 
    is_inner::Bool= false,
    do_hessian::Bool = false,
    do_guess_hist::Bool=false,
    lagrange::Bool = false,
    hesse_extra::Bool = false,
    basis_size::Int64=100,
    init_position::Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}}}=nothing)
    #this is a wrapper function to call the inner_procedure function with the correct arguments 
    # done such that the sequentiell code has a proper name
    return inner_procedure!(q_num=q_num,num_steps=num_steps, delta_gam=delta_gam,quiet=quiet,coarse_cacher=coarse_cacher,dist_break=dist_break,
    eps=eps,max_num_steps=max_num_steps,minier=minier,potential=potential,siO_glass=siO_glass,q_guesser_cache=q_guesser_cache,lagrange=lagrange, hesse_extra=hesse_extra, basis_size=basis_size,
    mini_breaker_search=mini_breaker_search,is_inner=is_inner,do_hessian=do_hessian,init_position=init_position, do_guess_hist=do_guess_hist)
end

function inner_procedure!(;
    q_num::Int64=0,
    num_steps::Int64= 35,
    delta_gam::Float64 = 0.005,
    quiet::Bool= false,
    coarse_cacher::Bool= false,
    dist_break::Union{Nothing,Float64}= nothing,
    eps::Float64= 1.0e-5,
    max_num_steps::Int64= 100000,
    minier::Union{Nothing,Minimiser}= nothing,
    potential::Union{Nothing,Int64}= nothing,
    siO_glass::Bool= false,
    q_guesser_cache::Union{Nothing,Vector{Matrix{Float64}}}= nothing,
    mini_breaker_search::Union{Nothing,Float64}= nothing, 
    is_inner::Bool= false,
    do_hessian::Bool = false,
    do_guess_hist::Bool=false,
    lagrange::Bool=false,
    hesse_extra::Bool=false,
    basis_size::Int64=100, 
    init_position::Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}}}=nothing)::Union{Nothing, Dict{String,Union{Float64, Minimiser, Vector{Float64}, Vector{Int64},Matrix{Float64}, Vector{Tuple{Vector{Float64},Vector{Float64}}},Vector{Matrix{Float64}}, Nothing}}}
    
    
    if hesse_extra
        lagrange = false
        q_num = 2
        reduced_basis_size = basis_size
        alpha_hesse = zeros(Float64, reduced_basis_size-2)
        molhess_object = Jumol.Hessian(minier.calc)
        Jumol.calc_hessian_anal!(molhess_object)
        eigens = eigen(molhess_object.hessian_mat)
        pairs_eig = Tuple{Float64, Vector{Float64}}[(eigens.values[j], eigens.vectors[:,j]) for j in 1:length(eigens.values)]
        pairs_eig = sort(pairs_eig, by=x->x[1])
        hesse_basis = Jumol.build_basis(pairs_eig, reduced_basis_size)
    end

    if isnothing(minier)
        molstruc_init, molcalc = create_struct(network_glass=siO_glass)
        Jumol.initialize_potential!(molcalc, potential)
        minier = Jumol.Minimiser(molcalc)
    else 
        if !isnothing(init_position)
            set_coordinates!(minier.calc.structure, init_position)
        end
    end
    
    if coarse_cacher
        coarse_cache = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, num_steps+1)
        Jumol.get_global_atom_pos_vec(minier.calc.structure)
        coarse_cache[1] = (Jumol.get_global_atom_pos_vec(minier.calc.structure),Jumol.get_box_config(minier.calc.structure))
    else
        coarse_cache = nothing
    end
    if isnothing(dist_break)
        limiter = 999.
    else
        limiter = dist_break - 5.0e-3
    end
    
    minier.acc_factor = 1.001
    minier.alpha_min = 1.0e0
    q_counter = 0
    if !isnothing(q_guesser_cache)
        cache_size = size(q_guesser_cache[1], 2)
        if cache_size != q_num
            q_guesser_cache = nothing
        else
            q_counter = q_num
        end
        if minier.cache_used
            init_guess = minier.init_guess
            after_p = minier.after_p
        else
            init_guess = Vector{Float64}(undef, 3*minier.calc.structure.noa) 
            after_p = Vector{Float64}(undef, 3*minier.calc.structure.noa)
        end
    elseif q_num > 1
        if minier.cache_used
            init_guess = minier.init_guess
            after_p = minier.after_p
        else
            init_guess = Vector{Float64}(undef, 3*minier.calc.structure.noa) 
            after_p = Vector{Float64}(undef, 3*minier.calc.structure.noa)
        end
    else
        init_guess = nothing
    end
    if minier.cache_used
        U_hist = minier.U_hist
        r_full = minier.r_full
        rt_full = minier.rt_full
        step_hist = minier.step_hist
        stress_hist = minier.stress_hist
        positions = minier.positions
    else
        U_hist = Vector{Float64}(undef, num_steps) 
        r_full = zeros(Float64, 3*minier.calc.structure.noa, num_steps) 
        rt_full = zeros(Float64, 3*minier.calc.structure.noa, num_steps) 
        step_hist = Vector{Int64}(undef, num_steps)
        stress_hist = Vector{Float64}(undef, num_steps)
        positions = Vector{Float64}(undef, 3*minier.calc.structure.noa)
        if do_guess_hist
            guess_hist = Matrix{Float64}(undef, 3*minier.calc.structure.noa, num_steps) 
        else
            guess_hist = nothing
        end
    end
    end_procedure = false   
    starting_time = time() ##
    total_time = time()
    if do_hessian
        hessian_hist = Vector{Matrix{Float64}}(undef, num_steps)
    else
        hessian_hist = nothing
    end
    dist_guess = Float64[]
    time_guessing = Float64[]
    final_dist = Float64[]
    for i in 1:num_steps
        if !(is_inner&& i==1) # first inner step does no affine displacement but only minimization
            Jumol.set_affine_deform_shear!(minier.calc.structure, delta_gam, 0.0, 0.0)
            #Jumol.set_affine_deform_vol!(Deformer, delta_gam, -0.0/10.0, 0.)
            Jumol.update_distances!(minier.calc.structure)
        end
        Jumol.get_positions_vector(minier.calc.structure, positions)
        r_full[:, i] .= positions
        if !isnothing(q_guesser_cache) && i < q_num + 1
            Jumol.push_matrix!(q_guesser_cache[1], positions)
        end
        if q_counter >= q_num && q_num > 1
            before_guess = time()
            if lagrange
                @views Jumol.lagrange_extrapolation!(init_guess, r_full[:,i-q_num+1:i], rt_full[:,i-q_num:i-1], i*delta_gam, delta_gam)
            elseif hesse_extra
                #println("Dist between old init and actual is ", norm(abs.(r_full[:,i-1] + hesse_basis*alpha_hesse - rt_full[:,i-1])), " and guess ", norm(abs.(init_guess - rt_full[:,i-1])), " a priori ", norm(abs.(init_guess - r_full[:,i-1])), " and real ", norm(abs.(rt_full[:,i-1] - r_full[:,i-1])))
                #println("Dist between old and new alpga is ", norm(abs.(alpha_hesse - best_approx_full(r_full[:,i-1], rt_full[:,i-1], hesse_basis))))
                #pretty_table(alpha_hesse)
                #pretty_table(best_approx_full(r_full[:,i-1], rt_full[:,i-1], hesse_basis))
                alpha_hesse .= best_approx_full(r_full[:,i-1], rt_full[:,i-1], hesse_basis)
                @views status = Jumol.hesse_extrapolation!(init_guess, r_full[:,i-1:i], hesse_basis, alpha_hesse)
                if !status
                    println("There is a platic event here")
                    molhess_object = Jumol.Hessian(minier.calc)
                    Jumol.calc_hessian_anal!(molhess_object)
                    eigens = eigen(molhess_object.hessian_mat)
                    pairs_eig = Tuple{Float64, Vector{Float64}}[(eigens.values[j], eigens.vectors[:,j]) for j in 1:length(eigens.values)]
                    pairs_eig = sort(pairs_eig, by=x->x[1])
                    hesse_basis .= Jumol.build_basis(pairs_eig, reduced_basis_size)
                end
            else
                if !isnothing(q_guesser_cache) && i < q_num + 1 # given a q cache from another processor we use that until our r_full matrix filled up enough
                    Jumol.construct_init_guess!(q_guesser_cache[1], q_guesser_cache[2], init_guess)
                else
                    @views Jumol.construct_init_guess!(r_full[:,i-q_num+1:i], rt_full[:, i-q_num:i-1], init_guess)
                end
            end
            push!(dist_guess, norm(abs.(init_guess - r_full[:,i])))
            push!(time_guessing, time()-before_guess)
            if do_guess_hist
                guess_hist[:,i] .= init_guess
            end
            #println("Dist between steps is ", norm(abs.(positions - init_guess)))
            Jumol.run_cg!(minier,num_steps = max_num_steps, tolerance = eps, init_guess = init_guess, quiet = true, dist_break = dist_break, final_stress = true, breaker_on_dist = mini_breaker_search)
        else
            Jumol.run_cg!(minier,num_steps = max_num_steps, tolerance = eps, init_guess = nothing, quiet = true, dist_break = dist_break, final_stress = true, breaker_on_dist = mini_breaker_search)
        end
        q_counter += 1
        temp_flag = -1.
        if minier.break_criterion=="Distance between steps"
            println("AQS step is breaking in step ", i)
        end
        if minier.alpha_min < 1.0e-4
            temp_flag = 99.
        end
        stress_hist[i] = sum(minier.calc.stress_tensor)
        if minier.is_broken
            end_procedure = true
        end

        step_hist[i] = minier.num_steps_real
        if q_num >1
            Jumol.get_positions_vector(minier.calc.structure, after_p)
            if !isnothing(q_guesser_cache) && i < q_num + 1 
                Jumol.push_matrix!(q_guesser_cache[2], after_p)
            end
        end

        if !isnothing(init_guess)
            if isnothing(mini_breaker_search)
                @views avg_dist_aff_nonaff = maximum([norm(abs.(init_guess - after_p)), temp_flag])
            else
                @views avg_dist_aff_nonaff = norm(abs.(positions - after_p))
            end

            if !quiet
                println("Absolute distance start guess to solution is ", avg_dist_aff_nonaff)
            end

            if avg_dist_aff_nonaff > limiter && isnothing(mini_breaker_search)
                if !quiet
                    #println("Distance to large doing it again")
                end
                Jumol.run_cg!(minier, num_steps = max_num_steps, tolerance = eps, init_guess = positions, quiet = true, final_stress = true)#, dist_break = dist_break)
                stress_hist[i] = sum(minier.calc.stress_tensor)
                Jumol.get_positions_vector(minier.calc.structure, after_p)
                step_hist[i] += minier.num_steps_real
                avg_dist_aff_nonaff = norm(abs.(positions - after_p))
                q_counter = 0
            end
        end
        
        Jumol.put_atoms_back_to_box!(minier.calc.structure) 
        ##Need to be added back later but this makes it harder to measure actuall distance between two differnet steps 
        if q_num > 1
            rt_full[:,i] .= after_p
        else
            Jumol.get_positions_vector(minier.calc.structure, positions)
            rt_full[:, i] .= positions
        end
        if q_counter >= q_num && q_num > 1 && !isnan(norm(abs.(init_guess - rt_full[:,i])))
            dd = norm(abs.(init_guess - rt_full[:,i]))
            if dd< 10.0
                push!(final_dist, dd)
            end
        end
        if !quiet
            println("U is ", minier.pot_en_hist[minier.num_steps_real])
        end
        U_hist[i] = minier.pot_en_hist[minier.num_steps_real]
        if coarse_cacher
            coarse_cache[i+1] = (Jumol.get_global_atom_pos_vec(minier.calc.structure), Jumol.get_box_config(minier.calc.structure))
        end
        if end_procedure
            break
        end
        if do_hessian
            molhess = Jumol.Hessian(minier.calc)
            Jumol.calc_hessian_anal!(molhess);
            hessian_hist[i] .= molhess.hessian_mat
        end
    end
    #println("Mean dist guess is ", mean(dist_guess))
    #println("Time spent guessing ", sum(time_guessing))
    #println("Mean final dist is ", mean(final_dist))
    stress_hist ./= minier.calc.structure.noa
    if minier.cache_used
        return nothing
    else
        return Dict("data"=>minier, "U_hist"=>U_hist, "step_hist"=>step_hist, "r_full"=>r_full, "rt_full"=>rt_full,
                "stress_hist"=>stress_hist, "hessian_hist"=>hessian_hist, "coarse_cache"=>coarse_cache, "guess_hist"=>guess_hist, "elapsed_time"=>time()-total_time)
    end
end

"""
    make_map_w_s(maxnum)
    Creates a map between which worker handles which remote parallel task.
"""
function make_map_w_s!(worker_step_map,maxnum::Int64, max_key::Int64)
    """
    creates map to know on which worker which level is running, allocation is random hence this is important to keep track where which solutions are held
    """;
    for o in 2:maxnum+1
        allow = @fetchfrom o @isdefined inner_marker
        if allow
            inner_marking = @fetchfrom o inner_marker
            if inner_marking <= max_key
                worker_step_map[inner_marking] = o
            end
        end
    end
end

"""
    extract_shared_array_to_dict(worker_mapping)
    Extracting from shared arrays all relevant information into an easy Dictonary. Better for readability and caching
"""
function extract_shared_array_to_dict!(info_dict:: Dict{Int64,Dict{String, Union{Array{Union{Float64, Int64}}, Nothing}}}, worker_mapping::Dict{Int64, Int64}; lazy_mode::Bool=false, lazy_index::Int64=-1)
    """
    the parallel calculation can only write into shared arrays and here we assemvel those to a dict again for easier readability and easier use of indices
    """;
    for c_step in eachindex(worker_mapping)
        if lazy_mode
            info_dict[c_step]["U"] = @fetchfrom worker_mapping[c_step] @views U
        else
            info_dict[c_step]["stress"] = @fetchfrom worker_mapping[c_step] @views stress_hist
            info_dict[c_step]["rs"] = @fetchfrom worker_mapping[c_step] @views r
            info_dict[c_step]["rts"] = @fetchfrom worker_mapping[c_step] @views rt
            info_dict[c_step]["steps"] = @fetchfrom worker_mapping[c_step] @views steps
        end
    end
end

function construct_dicts(max_coarse::Int64)
    info_dict = Dict{Int64,Dict{String, Union{Array{Union{Float64, Int64}}, Nothing}}}()
    for c_step in 1:max_coarse
        info_dict[c_step] = Dict{String, Union{Array{Union{Float64, Int64}}, Nothing}}()
        info_dict[c_step]["stress"] = nothing
        info_dict[c_step]["U"] = nothing
        info_dict[c_step]["rs"] = nothing
        info_dict[c_step]["rts"] = nothing
        info_dict[c_step]["steps"] = nothing
    end
    return info_dict, Dict{Int64, Int64}(k=>k for k in 1:max_coarse)
end

#########################################################
function aqs_sip_procedure(shifter::Aqs_parallel_shifter)
    shifter.track_progress = 0
    starter_all = time()
    time_of_parall = 0.
    U_tracker = Jumol.U_obj[]
    starter_all_list = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, shifter.P)
    info_dict, worker_step_map = construct_dicts(shifter.P)
    full_coordinate_history = Matrix{Float64}(undef, 3*shifter.minier.calc.structure.noa, shifter.aqs_steps)
    while shifter.track_progress < shifter.aqs_steps
        #println("Progress t=", shifter.track_progress, " time ", time()-starter_all, "s")
        coarse_range_left = Int64(ceil(minimum([shifter.P, (shifter.aqs_steps - shifter.track_progress)/shifter.D_split])))
        if coarse_range_left == 1
            last_steps = shifter.aqs_steps - shifter.track_progress
            start_last = time()
            dict_final_fine = inner_procedure!(q_num=shifter.q_num, num_steps=last_steps, 
                            eps= shifter.eps_mini, delta_gam=shifter.delta_gam, quiet=shifter.quiet, max_num_steps=shifter.max_min_steps, 
                            minier=shifter.minier, potential = shifter.potential, q_guesser_cache = shifter.q_cache)
            full_coordinate_history[:, shifter.track_progress+1:shifter.track_progress+last_steps] .= dict_final_fine["rt_full"][:,1:last_steps]
            
            time_of_parall += time()-start_last
            for final_acc in 1:last_steps
                @views u_temp = Jumol.U_obj(dict_final_fine["U_hist"][final_acc], norm(dict_final_fine["rt_full"][:, final_acc]), final_acc, 1, shifter.track_progress+final_acc,
                                    (1, shifter.track_progress), dict_final_fine["step_hist"][final_acc], dict_final_fine["stress_hist"][final_acc],  accepted = true)
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
                @sync @distributed for coarse_step in 1:coarse_range_left
                    #iii = @allocated begin
                        st = time()
                        if coarse_step == 1
                            inner_cache = shifter.q_cache #only makes sense to linear guess when we have fine steps before
                        else
                            inner_cache = nothing
                        end

                        parallel_steps = shifter.D_split
                        dict_fine =  Jumol.inner_procedure!(q_num=shifter.q_num, num_steps= parallel_steps + Int64(coarse_step!=1), coarse_cacher = false,
                                    eps= shifter.eps_mini, delta_gam=shifter.delta_gam, quiet=shifter.quiet, max_num_steps=shifter.max_min_steps, 
                                    minier=shifter.minier, potential=shifter.potential, q_guesser_cache=inner_cache, is_inner = coarse_step!=1,
                                    init_position=starter_all_list[coarse_step])
        
                        global parallel_steps
                        global inner_marker = coarse_step
                        global U = @views dict_fine["U_hist"]
                        global shared_position = (Jumol.get_global_atom_pos_vec(shifter.minier.calc.structure), Jumol.get_box_config(shifter.minier.calc.structure))
                        global stress_hist = dict_fine["stress_hist"]
                        global steps = dict_fine["step_hist"]
                        global r=dict_fine["r_full"]
                        global rt=dict_fine["rt_full"]
                    #end
                    #println("Alloc inside ", iii)
                end
            #end
            #println("The allocated in para was ", ddd)
            time_of_parall += time()-before_para
            #println("Para took ", time()-before_para)
            make_map_w_s!(worker_step_map, shifter.P, coarse_range_left) 
            extract_shared_array_to_dict!(info_dict, worker_step_map, lazy_mode = true);
            choosing_time = time()

            # choosing sequence starts here
            for check_index in 1:coarse_range_left-1
                if check_index==1
                    ending = shifter.D_split
                else
                    ending = shifter.D_split + 1
                end
                if abs(info_dict[check_index]["U"][ending] - info_dict[check_index+1]["U"][1]) > shifter.eps_accepter || shifter.track_progress + check_index*shifter.D_split >= shifter.aqs_steps
                    set_coordinates!(shifter.minier.calc.structure, @fetchfrom worker_step_map[check_index] shared_position)
                    extract_shared_array_to_dict!(info_dict, worker_step_map, lazy_mode = false);
                    shifter.q_cache = nothing #because linear guessing is reset when plastic event happens
                    Jumol.push_dict_to_U!(U_tracker, info_dict, shifter.track_progress, shifter.D_split, check_index)
                    update_rts_hist!(full_coordinate_history, shifter, info_dict, check_index)
                    shifter.track_progress += check_index * shifter.D_split
                    break
                elseif check_index+1 == coarse_range_left
                    set_coordinates!(shifter.minier.calc.structure, @fetchfrom worker_step_map[check_index+1] shared_position)
                    extract_shared_array_to_dict!(info_dict, worker_step_map, lazy_mode = false);
                    shifter.q_cache = @views Matrix{Float64}[info_dict[check_index+1]["rs"][:, shifter.D_split-shifter.q_num+1:shifter.D_split], info_dict[check_index+1]["rts"][:, shifter.D_split-shifter.q_num+1:shifter.D_split]]
                    Jumol.push_dict_to_U!(U_tracker, info_dict, shifter.track_progress, shifter.D_split, check_index+1)
                    update_rts_hist!(full_coordinate_history, shifter, info_dict, check_index)
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

function update_rts_hist!(rts_hist::Matrix{Float64}, shifter::Aqs_parallel_shifter, info_dict::Dict{Int64,Dict{String, Union{Array{Union{Float64, Int64}}, Nothing}}}, boarder::Int64)
    cur_t = shifter.track_progress+1
    #flawed since it does not take into account if there was a shift or not, which is flagged as is_inner
    for coarse_k in 1:maximum(keys(info_dict))
        if coarse_k<=boarder
            this_size = size(info_dict[coarse_k]["rts"], 2)
            if this_size > shifter.D_split
                this_size = shifter.D_split
            end
            @views rts_hist[:, cur_t:cur_t+this_size-1] .= info_dict[coarse_k]["rts"][:, end-this_size+1:end]
            cur_t += this_size
        end
    end
end

function set_coordinates!(structure::Structure, init_position::Tuple{Vector{Float64}, Vector{Float64}})
    set_global_pos_vec!(structure, init_position[1])
    set_box_config!(structure, init_position[2])
    set_box_basis_vectors!(structure.box)
    update_distances!(structure)
    if structure.linked_cells_bool == true
        update_cell!(structure)
        construct_neigh_cell_list!(structure)
        linked_cell_list!(structure)
    end
end