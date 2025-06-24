## minimizer
######################################
######################################
######################################


## struct Min definition
mutable struct Minimiser
    calc::Calc
    alpha_min::Float64
    full_update_num::Int64
    acc_factor::Float64
    output_filename::String
    num_steps_max::Int64
    num_steps_real::Int64
    is_broken::Bool
    break_criterion::String
    pot_en_hist::Vector{Float64}
    force_hist::Vector{Float64}
    contr_vec::Vector{Int64}
    current_pos::Vector{Float64}
    placeholder_force_vec::Vector{Float64}
    rp1::Vector{Float64}
    r::Vector{Float64}
    F::Vector{Float64}
    Fp1::Vector{Float64}
    p::Vector{Float64}
    delta_step_vec::Vector{Float64}
    group_list::Vector{Int64}
 
    cache_used::Bool
    D_split::Int64
    init_guess::Union{Nothing,Vector{Float64}}
    after_p::Union{Nothing,Vector{Float64}}
    U_hist::Union{Nothing,Vector{Float64}}
    r_full::Union{Nothing,Matrix{Float64}}
    rt_full::Union{Nothing,Matrix{Float64}}
    step_hist::Union{Nothing,Vector{Int64}}
    stress_hist::Union{Nothing,Vector{Float64}}
    positions::Union{Nothing,Vector{Float64}}
    
    function Minimiser(calc;
        alpha_min=1.0e-2,
        full_update_num=10,
        acc_factor=1.01,
        output_filename="None",
        num_steps_max = 100000,
        num_steps_real = -1,
        is_broken = false,
        break_criterion = "No termination",
        pot_en_hist = zeros(Float64, num_steps_max),
        force_hist = zeros(Float64, num_steps_max),
        contr_vec = zeros(Int64, 3*calc.structure.noa),
        current_pos = Vector{Float64}(undef,3*calc.structure.noa),
        placeholder_force_vec = Vector{Float64}(undef,3*calc.structure.noa),
        rp1 = Vector{Float64}(undef,3*calc.structure.noa),
        r = Vector{Float64}(undef,3*calc.structure.noa),
        F = Vector{Float64}(undef,3*calc.structure.noa),
        Fp1 = Vector{Float64}(undef,3*calc.structure.noa),
        p = Vector{Float64}(undef,3*calc.structure.noa),
        delta_step_vec = Vector{Float64}(undef, 3*calc.structure.noa),
        group_list = Int64[0],
        cache_used=false, D_split = 0) 
        if cache_used
            init_guess = Vector{Float64}(undef, 3*calc.structure.noa) 
            after_p = Vector{Float64}(undef, 3*calc.structure.noa) 
            U_hist = Vector{Float64}(undef, D_split+1)
            r_full = Matrix{Float64}(undef, 3*calc.structure.noa,D_split+1)
            rt_full = Matrix{Float64}(undef, 3*calc.structure.noa,D_split+1)
            step_hist =  Vector{Int64}(undef, D_split+1)
            stress_hist =  Vector{Float64}(undef, D_split+1)
            positions = Vector{Float64}(undef, 3*calc.structure.noa)
        else
            init_guess = nothing
            after_p = nothing
            U_hist = nothing
            r_full = nothing
            rt_full = nothing
            step_hist = nothing
            stress_hist = nothing
            positions = nothing
        end
        return new(calc,alpha_min,full_update_num,acc_factor,output_filename,num_steps_max,
                    num_steps_real,is_broken,break_criterion, pot_en_hist, force_hist, contr_vec, 
                    current_pos,placeholder_force_vec, rp1, r, F, Fp1, p, delta_step_vec, group_list,
                    cache_used, D_split, init_guess, after_p, U_hist, r_full, rt_full, step_hist, stress_hist, positions)
    end
end

## run steepest decent algorithm
function run_sd!(minier::Minimiser; num_steps::Int64=100000, tolerance =1.0e-7, init_guess=nothing, dist_break = nothing, full_update=true, save_lammpstrj=false, save_step=100, quiet = true)
    ##
    #if save_lammpstrj == true
     #   open_file!(minier.calc.structure.reader, minier.output_filename, "w")
    #end
    if minier.num_steps_max != num_steps
        minier.num_steps_max = num_steps
        minier.pot_en_hist = zeros(Float64, num_steps)
        minier.force_hist = zeros(Float64, num_steps)
    else
        minier.pot_en_hist .= 0.0
        minier.force_hist .= 0.0
    end
    ##
    max_disp_step = 0.0
    ## full update before the minimization starts
    if full_update==true
        if minier.calc.structure.linked_cells_bool == true
            update_cell!(minier.calc.structure) # needed updating as structure was not defined 
        	linked_cell_list!(minier.calc.structure)
        	construct_neigh_cell_list!(minier.calc.structure)
        end
        update_distances!(minier.calc.structure)
    end
    ## control vector for tolerance
    if minier.calc.structure.all_grous_zero
        minier.contr_vec .= 1
    else
        minier.contr_vec .= 0
        cntr = 1
        for i in 1:minier.calc.structure.noa
            for ii in eachindex(minier.group_list)
                if minier.calc.structure.atom_list[i].group == minier.group_list[ii]
                    minier.contr_vec[cntr:cntr+2] .= 1
                end
            end
            cntr += 3
        end
    end

    if !isnothing(init_guess)
        for i in eachindex(minier.calc.structure.atom_list)
            minier.calc.structure.atom_list[i].pos = init_guess[(i-1)*3+1:i*3]
        end
    end

    U = calc_total_pot_en!(minier.calc)
    placeholder_pos = @MVector zeros(Float64, 3)
    placeholder_t = @MVector zeros(Float64, 3)
    placeholder_force = @MVector zeros(Float64, 3)
    minier.placeholder_force_vec = zeros(Float64, 3*minier.calc.structure.noa)
    ## run through the number of minimization steps
    for i in 1:num_steps
        ## force calculation
        calc_all_pair_forces!(minier.calc)
        ## add the step downwards to the position of every atom
        max_delta_step = 0.0
        current_pos = zeros(Float64, 3*minier.calc.structure.noa)

        for ii in 1:minier.calc.structure.noa
            current_pos[(ii-1)*3+1:ii*3] = copy(minier.calc.structure.atom_list[ii].pos)
            for iii in eachindex(minier.group_list)
                if minier.calc.structure.atom_list[ii].group == minier.group_list[iii]
                    delta_step = minier.calc.structure.atom_list[ii].force*minier.alpha_min
                    minier.calc.structure.atom_list[ii].pos += delta_step
                    len_step = norm(delta_step)
                    if len_step > max_delta_step
                        max_delta_step = len_step
                    end
                end
            end
        end

        if !isnothing(init_guess) && !isnothing(dist_break)
            disti = norm(init_guess - current_pos) 
            #println(disti, dist_break)
            if disti > dist_break 
                minier.num_steps_real = i
                #println("Breaking out because to close, steps ", i)
                @goto escape_loop
            end
        end
        ## update distances
        max_disp_step += max_delta_step
        if max_disp_step > minier.calc.structure.rskin*0.5
            if full_update==true
                if minier.calc.structure.linked_cells_bool == true
                    update_cell!(minier.calc.structure)
                	linked_cell_list!(minier.calc.structure)
                	construct_neigh_cell_list!(minier.calc.structure)
                end
                update_distances!(minier.calc.structure)
                max_disp_step = 0.0
            else
                update_distances_partly!(minier.calc.structure)
            end
        else
            update_distances_partly!(minier.calc.structure)
        end
        ## potential calculation
        U = calc_total_pot_en!(minier.calc)
        if save_lammpstrj == true
            if (i-1)%save_step == 0
                write_box!(minier.calc.structure.reader,minier.calc.structure,i)
            end
        end
        minier.pot_en_hist[i] = U
        ## force_tolerance achieved -> break
        F = get_global_force_vec(minier.calc.structure, minier.placeholder_force_vec) .* minier.contr_vec
        f_scal = dot(F,F)
        minier.force_hist[i] = f_scal
        if !quiet
            print("\e[2K")
            print("\r potential: ", U, " force norm: ", f_scal, " alpha_min: ", minier.alpha_min)
            println("This is the f_scal for term or not", f_scal)
        end
        if f_scal < tolerance
            if !quiet
                println("\nforce tolerance achieved at step: ", i)
            end 
            minier.num_steps_real = i
            break
        end
        ## update acutal number of steps done
        minier.num_steps_real = i
    end
    @label escape_loop
    if save_lammpstrj == true
        close_file!(minier.calc.structure.reader)
    end
end

function shift_correct(incoming::Vector{Float64}; r_full::Union{Vector{Float64},Nothing}=nothing)
    if isnothing(r_full)
        min_end = 0.0
    else
        min_end = minimum([r_full[k] for k in eachindex(r_full) if mod(k,3) != 0])
    end

    outgoing = copy(incoming) .- minimum(incoming) .+ min_end
    size_d = Int64(size(incoming,1)/3)
    for i in 1:size_d
        outgoing[3*i] = 0.0
    end
    return outgoing
end

function run_cg!(minier::Minimiser;
                num_steps::Int64 = 100000,
                tolerance::Float64 = 1.0e-4,
                init_guess::Union{Nothing, Vector{Float64}} = nothing,
                dist_break::Union{Nothing, Float64} = nothing,
                full_update::Bool = true,
                save_lammpstrj::Bool = false,
                save_step::Int64 = 100,
                quiet::Bool = true,
                final_stress::Bool = false,
                breaker_on_dist::Union{Nothing, Float64} = nothing)

    if num_steps==0
        minier.num_steps_real = 1
        return nothing
    end
    ##
    if save_lammpstrj == true
        open_file!(minier.calc.structure.reader, minier.output_filename, "w")
    end
    ##
    minier.is_broken = false
    if minier.num_steps_max != num_steps
        minier.num_steps_max = num_steps
        minier.pot_en_hist = zeros(Float64, num_steps)
        minier.force_hist = zeros(Float64, num_steps)
    else
        minier.pot_en_hist .= 0.0
        minier.force_hist .= 0.0
    end
    ## full update before the minimization starts
    max_disp_step = 0.0
    #println("Minier pot ", minier.calc.potential)
    #println("Minier cur pot ", get_current_pot(minier.calc.structure))
    if !isnothing(init_guess)
        for k in 1:minier.calc.structure.noa
            @views minier.calc.structure.atom_list[k].pos .= init_guess[3*(k-1)+1 : 3*k]
        end 
    end
    
    if full_update
        if minier.calc.structure.linked_cells_bool == true
            update_cell!(minier.calc.structure)
            linked_cell_list!(minier.calc.structure)
            construct_neigh_cell_list!(minier.calc.structure)
        end
        update_distances!(minier.calc.structure)
    end
    ## minimization contribution vector (minimize only those atoms that are part of a certain group)
    if minier.calc.structure.all_grous_zero
        minier.contr_vec .= 1
    else
        minier.contr_vec .= 0
        cntr = 1
        for i in 1:minier.calc.structure.noa
            for ii in eachindex(minier.group_list)
                if minier.calc.structure.atom_list[i].group == minier.group_list[ii]
                    minier.contr_vec[cntr:cntr+2] .= 1
                end
            end
            cntr += 3
        end
    end
    c1 = 0.0001
    c2 = 0.1
    ## initial calculation
    Up1 = 0.0 
    get_global_atom_pos_vec!(minier.calc.structure, minier.r)
    U = calc_total_pot_en!(minier.calc)
    init_U = U 
    calc_all_pair_forces!(minier.calc)
    minier.F .= get_global_force_vec(minier.calc.structure, minier.F)
    minier.F .*= minier.contr_vec
    minier.Fp1 .= 0.0
    minier.p .= minier.F
    cg_temp_alpha = minier.alpha_min
    minier.current_pos .= 0.0
    d2_cache_valid = false
    d2_cache = 0.0
    ## run minimization loop
    for i in 1:num_steps
        ## testing alpha using the WOLFE conditions
        WOLFE_TOTAL = false
        minier.current_pos .= 0.0
        while WOLFE_TOTAL == false
            WOLFE1 = false
            WOLFE2 = false
            minier.delta_step_vec .= cg_temp_alpha.*minier.p
            minier.rp1 .= minier.r .+ minier.delta_step_vec
            ## update positions, calculate potential engergy and forces
            set_global_pos_vec!(minier.calc.structure, minier.rp1)
            put_atoms_back_to_box!(minier.calc.structure)

            if max_disp_step > minier.calc.structure.rskin*0.5 # half of the skin distance because center atom can move in oposite direction
                if full_update==true
                    if minier.calc.structure.linked_cells_bool == true
                        update_cell!(minier.calc.structure)
                    	linked_cell_list!(minier.calc.structure)
                    	construct_neigh_cell_list!(minier.calc.structure)
                    end
                    update_distances!(minier.calc.structure)
                    max_disp_step = 0.0
                else
                    update_distances_partly!(minier.calc.structure)
                end
            else
                update_distances_partly!(minier.calc.structure)
            end
            Up1 = calc_total_pot_en!(minier.calc)
            calc_all_pair_forces!(minier.calc)
            @views minier.Fp1 .= get_global_force_vec(minier.calc.structure, minier.placeholder_force_vec) .* minier.contr_vec
            if !d2_cache_valid
                d2_cache = dot(minier.p,minier.F)
                d2_cache_valid = true
            end
            @views checker = U - c1*cg_temp_alpha*d2_cache

            if Up1 <= checker
                WOLFE1 = true
            else
                WOLFE1 = false
            end
           
            if WOLFE1
                d1 = dot(minier.Fp1,minier.p)
                if d1 > c2*d2_cache
                    WOLFE2 = true
                else
                    WOLFE2 = false
                end
            end
            minier.pot_en_hist[i] = Up1 #corrected mistake, always save potential energy no matter if wolfe condiitions are fullfilled or not 
            ## check possibilities and adapt alpha_min accordingly
            if WOLFE1 == true && WOLFE2 == true
                #println("Both wolfes true")
                WOLFE_TOTAL = true
                d2_cache_valid = false
                ## Polak-Ribiere for the next step
                d3 = dot(minier.Fp1, minier.Fp1)
                @views beta = (d3 - dot(minier.Fp1, minier.F)) / dot(minier.F,minier.F) # sollte das nicht dot(minier.Fp1,F) sein ?
                @views minier.p .*= beta
                @views minier.p .+= minier.Fp1
                d2_cache_valid = true
                @views d2_cache = d3 + beta*d1
                minier.r .= minier.rp1
                U = Up1
                minier.F .= minier.Fp1
                cg_temp_alpha = cg_temp_alpha*minier.acc_factor
                ## check maximum distance of the step and add it to the value
                max_disp = 0.0
                for ii in 1:minier.calc.structure.noa
                    minier.current_pos[(ii-1)*3+1:ii*3] .= minier.calc.structure.atom_list[ii].pos
                    max_disp_check = ( minier.delta_step_vec[3*ii-2]*minier.delta_step_vec[3*ii-2] +
                                    minier.delta_step_vec[3*ii-1]*minier.delta_step_vec[3*ii-1] +
                                    minier.delta_step_vec[3*ii-0]*minier.delta_step_vec[3*ii-0]  )
                    if max_disp < max_disp_check
                        max_disp = max_disp_check
                    end
                end

                if !isnothing(init_guess) && !isnothing(dist_break)
                    @views disti = norm(minier.current_pos-init_guess)
                    if dist_break < disti < minier.calc.structure.box.h1[1] #box length
                        minier.num_steps_real = i
                        println("Breaking out because to close, steps ", i, " with dist ", disti)
                        minier.break_criterion = "Distance between steps"
                        @goto escape_loop_cg
                    end
                end

                max_disp_step += sqrt(max_disp)
                ## save the potential energy
                #minier.pot_en_hist[i] = U
                ## save the scalar force
                minier.force_hist[i] = d3 #cache is acutally dot(Fp, Fp) above F update in line 391
                ## save lammpstrj if necessary
                if save_lammpstrj == true
                    if (i-1)%save_step == 0
                        write_box!(minier.calc.structure.reader,minier.calc.structure,i)
                    end
                end
            else
                WOLFE_TOTAL = false
                #println("step: ", i)
                #println("WOLFE1: ", WOLFE1)
                #println("WOLFE2: ", WOLFE2)
                ## decrease step size
                cg_temp_alpha = cg_temp_alpha*0.5
                #println("alpha_min: ", cg_temp_alpha)
                ## search alpha too small -> break
                if cg_temp_alpha < 1e-4
                    if !quiet
                        println(" ")
                        println("break criterion alpha_min (inside): ", i," pot ", Up1)
                    end
                    minier.break_criterion = "Alpha min tolerance reached 1"
                    ## update acutal number of steps done
                    #println("Broken because of alpha_min")
                    minier.num_steps_real = i
                    break
                end
            end
        end
        ## search alpha too small -> break
        if cg_temp_alpha < 1e-4
            if !quiet
                println(" ")
                println("break criterion alpha_min: ", i," pot ", Up1)
            end
            ## update acutal number of steps done
            minier.num_steps_real = i
            minier.break_criterion = "Alpha min tolerance reached 2"
            #println("Broken because of alpha_min")
            break
        end

        ## force_tolerance achieved -> break
        force_norm = dot(minier.F,minier.F)
        if !quiet
            print("\e[2K \e[2K")
            print("\r potential: ", U, " force norm: ", force_norm, " alpha_min: ", cg_temp_alpha)
        end

        if force_norm < tolerance
            if !quiet
                println(" ")
                println("\nforce tolerance achieved at step: ", i)
            end
            #println("Broken because of force norm ", force_norm)
            ## update acutal number of steps done
            minier.num_steps_real = i
            minier.break_criterion = "Force norm < tol"
            break
        end
        ## update acutal number of steps done
        minier.num_steps_real = i

        if !isnothing(breaker_on_dist) && abs(init_U -  minier.pot_en_hist[i]) > breaker_on_dist
            #println("Breaking minimisation in coarse solver")
            minier.is_broken = true
            minier.break_criterion = "Break in coarse solver"
            break
        end
    end
    if minier.num_steps_real == num_steps
        minier.break_criterion = "Maximum number of steps reached"
    end
    @label escape_loop_cg
    if final_stress
        calc_only_stress!(minier.calc)
    end
    if save_lammpstrj == true
        close_file!(minier.calc.structure.reader)
    end
    return nothing
end

## run cg in reduced space
function run_cg_red_space!(minier::Minimiser, r_0::Vector{Float64}, space::Matrix{Float64};
                num_steps::Int64 = 100000,
                tolerance::Float64 = 1.0e-4,
                full_update::Bool = true,
                save_step::Int64 = 100,
                quiet::Bool = true,
                final_stress::Bool = false)

    minier.is_broken = false
    if minier.num_steps_max != num_steps
        minier.num_steps_max = num_steps
        minier.pot_en_hist = zeros(Float64, num_steps)
        minier.force_hist = zeros(Float64, num_steps)
    else
        minier.pot_en_hist .= 0.0
        minier.force_hist .= 0.0
    end
    ## full update before the minimization starts
    max_disp_step = 0.0
    
    if full_update
        if minier.calc.structure.linked_cells_bool == true
            update_cell!(minier.calc.structure)
            linked_cell_list!(minier.calc.structure)
            construct_neigh_cell_list!(minier.calc.structure)
        end
        update_distances!(minier.calc.structure)
    end
    ## minimization contribution vector (minimize only those atoms that are part of a certain group)
    if minier.calc.structure.all_grous_zero
        minier.contr_vec .= 1
    else
        minier.contr_vec .= 0
        cntr = 1
        for i in 1:minier.calc.structure.noa
            for ii in eachindex(minier.group_list)
                if minier.calc.structure.atom_list[i].group == minier.group_list[ii]
                    minier.contr_vec[cntr:cntr+2] .= 1
                end
            end
            cntr += 3
        end
    end
    c1 = 0.0001
    c2 = 0.1
    ## initial calculation
    red_r = zeros(Float64, size(space,2))
    red_rp1 = zeros(Float64, size(space,2))
    red_step_vec = zeros(Float64, size(space,2))
    z_vec = zeros(Float64, 3*minier.calc.structure.noa)
    Up1 = 0.0 
    U = calc_total_pot_en!(minier.calc)
    init_U = U 
    calc_all_pair_forces!(minier.calc)
    
    red_F = transpose(space) * get_global_force_vec(minier.calc.structure, minier.F) 
    red_Fp1 = zeros(Float64, size(space,2))
    red_p = red_F
    cg_temp_alpha = minier.alpha_min
    minier.current_pos .= 0.0
    d2_cache_valid = false
    d2_cache = 0.0
    ## run minimization loop
    for i in 1:num_steps
        ## testing alpha using the WOLFE conditions
        WOLFE_TOTAL = false
        minier.current_pos .= 0.0
        while WOLFE_TOTAL == false
            WOLFE1 = false
            WOLFE2 = false
            red_step_vec .= cg_temp_alpha .* red_p
            red_rp1 .= red_r .+ red_step_vec
            z_vec .= r_0 .+ space * red_rp1
            ## update positions, calculate potential engergy and forces
            set_global_pos_vec!(minier.calc.structure, z_vec)
            put_atoms_back_to_box!(minier.calc.structure)

            if max_disp_step > minier.calc.structure.rskin*0.5 # half of the skin distance because center atom can move in oposite direction
                if full_update==true
                    if minier.calc.structure.linked_cells_bool == true
                        update_cell!(minier.calc.structure)
                    	linked_cell_list!(minier.calc.structure)
                    	construct_neigh_cell_list!(minier.calc.structure)
                    end
                    update_distances!(minier.calc.structure)
                    max_disp_step = 0.0
                else
                    update_distances_partly!(minier.calc.structure)
                end
            else
                update_distances_partly!(minier.calc.structure)
            end
            Up1 = calc_total_pot_en!(minier.calc)
            calc_all_pair_forces!(minier.calc)
            red_Fp1 .= transpose(space) * get_global_force_vec(minier.calc.structure, minier.placeholder_force_vec) .* minier.contr_vec
            if !d2_cache_valid
                d2_cache = dot(red_p, red_F)
                d2_cache_valid = true
            end
            @views checker = U - c1*cg_temp_alpha*d2_cache

            if Up1 <= checker
                WOLFE1 = true
            else
                WOLFE1 = false
            end
           
            if WOLFE1
                d1 = dot(red_Fp1,red_p)
                if d1 > c2*d2_cache
                    WOLFE2 = true
                else
                    WOLFE2 = false
                end
            end
            minier.pot_en_hist[i] = Up1 #corrected mistake, always save potential energy no matter if wolfe condiitions are fullfilled or not 
            ## check possibilities and adapt alpha_min accordingly
            if WOLFE1 == true && WOLFE2 == true
                #println("Both wolfes true")
                WOLFE_TOTAL = true
                d2_cache_valid = false
                ## Polak-Ribiere for the next step
                d3 = dot(red_Fp1, red_Fp1)
                @views beta = (d3 - dot(red_Fp1, red_F)) / dot(red_F,red_F)
                @views red_p .*= beta
                @views red_p .+= red_Fp1
                d2_cache_valid = true
                @views d2_cache = d3 + beta*d1
                red_r .= red_rp1
                U = Up1
                red_F .= red_Fp1
                cg_temp_alpha = cg_temp_alpha*minier.acc_factor
                ## check maximum distance of the step and add it to the value
                max_disp = 0.0
                for ii in 1:minier.calc.structure.noa
                    minier.current_pos[(ii-1)*3+1:ii*3] .= minier.calc.structure.atom_list[ii].pos
                    max_disp_check = ( minier.delta_step_vec[3*ii-2]*minier.delta_step_vec[3*ii-2] +
                                    minier.delta_step_vec[3*ii-1]*minier.delta_step_vec[3*ii-1] +
                                    minier.delta_step_vec[3*ii-0]*minier.delta_step_vec[3*ii-0]  )
                    if max_disp < max_disp_check
                        max_disp = max_disp_check
                    end
                end

                ## save the potential energy
                #minier.pot_en_hist[i] = U
                ## save the scalar force
                minier.force_hist[i] = d3 #cache is acutally dot(Fp, Fp) above F update in line 391
                ## save lammpstrj if necessary
                if save_lammpstrj == true
                    if (i-1)%save_step == 0
                        write_box!(minier.calc.structure.reader,minier.calc.structure,i)
                    end
                end
            else
                WOLFE_TOTAL = false
                #println("step: ", i)
                #println("WOLFE1: ", WOLFE1)
                #println("WOLFE2: ", WOLFE2)
                ## decrease step size
                cg_temp_alpha = cg_temp_alpha*0.5
                #println("alpha_min: ", cg_temp_alpha)
                ## search alpha too small -> break
                if cg_temp_alpha < 1e-4
                    if !quiet
                        println(" ")
                        println("break criterion alpha_min (inside): ", i," pot ", Up1)
                    end
                    minier.break_criterion = "Alpha min tolerance reached 1"
                    ## update acutal number of steps done
                    #println("Broken because of alpha_min")
                    minier.num_steps_real = i
                    break
                end
            end
        end
        ## search alpha too small -> break
        if cg_temp_alpha < 1e-4
            if !quiet
                println(" ")
                println("break criterion alpha_min: ", i," pot ", Up1)
            end
            ## update acutal number of steps done
            minier.num_steps_real = i
            minier.break_criterion = "Alpha min tolerance reached 2"
            #println("Broken because of alpha_min")
            break
        end

        ## force_tolerance achieved -> break
        force_norm = dot(red_F, red_F) #TODO
        if !quiet
            print("\e[2K \e[2K")
            print("\r potential: ", U, " force norm: ", force_norm, " alpha_min: ", cg_temp_alpha)
        end

        if force_norm < tolerance
            if !quiet
                println(" ")
                println("\nforce tolerance achieved at step: ", i)
            end
            #println("Broken because of force norm ", force_norm)
            ## update acutal number of steps done
            minier.num_steps_real = i
            minier.break_criterion = "Force norm < tol"
            break
        end
        ## update acutal number of steps done
        minier.num_steps_real = i
    end
    if minier.num_steps_real == num_steps
        minier.break_criterion = "Maximum number of steps reached"
    end
    @label escape_loop_cg
    if final_stress
        calc_only_stress!(minier.calc)
        #calc_stress_classic!(minier.calc)
    end

    return red_rp1, z_vec
end


function inverter(nummer::Float64)
    if abs(nummer)<1.0e-7
        return 0.0
    else
        return 1.0/nummer
    end
end

function nvt_steps!(min::Minimiser)
    if minier.calc.structure.linked_cells_bool == true
        update_cell!(minier.calc.structure)
        linked_cell_list!(minier.calc.structure)
        construct_neigh_cell_list!(minier.calc.structure)
    end
    update_distances!(minier.calc.structure)
    #
    num_steps = 100
    molunits = Jumol.Units(unit_type = 0)
    Int = Integrator(minier.calc, molunits)
    run_nvt!(Int, num_steps, 0.5, 0.5, 2.)
end
