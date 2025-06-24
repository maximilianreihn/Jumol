## force and potential calculations
######################################
######################################
######################################

include("pot/Yukawa_2D.jl")
include("pot/Tersoff.jl")
include("pot/LennardJones_MF.jl")
include("pot/Harm_ring_ring.jl")

mutable struct Calc
    structure::Structure
    type_potential::Int64
    placeholder_force::MVector{3, Float64}
    Calc(structure; type_potential = 99, placeholder_force = MVector{3,Float64}([0.,0.,0.])) = new(structure, type_potential, placeholder_force)
    potential::Union{Harm_ring_ring, Yukawa_2D, Tersoff, LJ}
    stress_tensor::Matrix{Float64}
    type_matrix::SMatrix{Int64}
end

## initialize type of potential before the calculation
function initialize_potential!(calc::Calc, type_potential::Int64)
    calc.type_potential = type_potential
    if calc.type_potential == 7
        calc.potential = Harm_ring_ring()
        set_ring_ring_params!(calc.potential)
    elseif calc.type_potential == 13
        calc.potential = Yukawa_2D(calc.structure.rc)
        set_Yukawa_2D_params!(calc.potential)
    elseif calc.type_potential == 21
        calc.potential = Tersoff(calc.structure.rc)
        set_Tersoff_params!(calc.potential, calc.structure)
    elseif calc.type_potential == 33
        calc.potential = LJ(calc.structure.rc)
        set_LJ_params!(calc.potential)
    end
    ## initialize stress tensor
    calc.stress_tensor = zeros(Float64, 3,3)
end


## set forces on all atoms to zero
function set_atom_forces_to_zero!(calc::Calc)
    for a in calc.structure.atom_list
        a.force .= 0.0
    end
    ## set also stress tensor to zero
    calc.stress_tensor .= 0.0
end

## set potentials of all atoms to zero
function set_atom_potential_to_zero!(calc::Calc)
    for a in calc.structure.atom_list
        a.pot=0.0
    end
end

## calculate the sum of all pair potentials on all atoms
function calc_all_pair_pot!(calc::Calc)
    set_atom_potential_to_zero!(calc)
    for (i,ii) in @views calc.structure.index_cache[1:calc.structure.index_counter]
        pot = 0.0
        #calc_pot!(calc,i,ii)
        if calc.type_potential == 7 || calc.type_potential == 13
            pot = calc_pot(calc.potential,calc.structure.atom_list[i].type,
                        calc.structure.atom_list[ii].type,
                        calc.structure.atom_list[i].norms[ii-i])
            calc.structure.atom_list[i].pot += pot
            calc.structure.atom_list[ii].pot += pot
        elseif calc.type_potential == 33
            pot = calc_pot(calc.potential, calc.structure.atom_list[i].type,
                        calc.structure.atom_list[ii].type,
                        calc.structure.atom_list[i].norms[ii-i])
            calc.structure.atom_list[i].pot += pot
            calc.structure.atom_list[ii].pot += pot
        elseif calc.type_potential == 21
            pot = calc_pot(calc.potential,i,ii)
            calc.structure.atom_list[i].pot += pot
        end
    end
end

## calculate total potential energy of the system
function calc_total_pot_en!(calc::Calc)
    calc_all_pair_pot!(calc)
    total_pot = 0.0
    for ii in 1:calc.structure.noa
        total_pot += calc.structure.atom_list[ii].pot
    end
    return total_pot
end

function get_current_pot(structure::Structure)
    total_pot = 0.0
    for ii in 1:structure.noa
        total_pot += structure.atom_list[ii].pot
    end
    return total_pot
end

## calculate the sum of all forces on each atom
function calc_all_pair_forces!(calc::Calc; calc_stress_tensor::Bool=false)
    set_atom_forces_to_zero!(calc)
    for (i,ii) in @views calc.structure.index_cache[1:calc.structure.index_counter]
        calc.placeholder_force .= 0.0
        calc_force!(calc,i,ii)
        calc.structure.atom_list[i].force += calc.placeholder_force
        calc.structure.atom_list[ii].force -= calc.placeholder_force
        if calc_stress_tensor
            if calc.structure.atom_list[i].calc_stress
                @views calc.stress_tensor .+= @views calc.structure.atom_list[i].distances[:,ii-i]*calc.placeholder_force'
            end
        end
    end
    if calc_stress_tensor
        calc.stress_tensor = calc.stress_tensor*(-1.0)/(calc.structure.box.lx*calc.structure.box.ly*calc.structure.box.lz)
    end
end

function calc_only_stress!(calc::Calc)
    for (i,ii) in @views calc.structure.index_cache[1:calc.structure.index_counter]
        calc.placeholder_force .= 0.0
        calc_force!(calc,i,ii)
        if calc.structure.atom_list[i].calc_stress
            @views calc.stress_tensor .+= @views calc.structure.atom_list[i].distances[:, ii-i] .* calc.placeholder_force'  # In-place addition of stress tensor
        end
    end
    calc.stress_tensor .*= -1.0 / (calc.structure.box.lx * calc.structure.box.ly * calc.structure.box.lz)
end

function calc_stress_classic!(calc::Calc)
    calc.stress_tensor .= 0.0
    #calc_all_pair_forces!(calc)
    for i in 1:calc.structure.noa
            @views calc.stress_tensor .+= calc.structure.atom_list[i].pos .* calc.structure.atom_list[i].force'
    end
    calc.stress_tensor .*= -1.0 / (calc.structure.box.lx * calc.structure.box.ly * calc.structure.box.lz)
end

## calculate the force for Tersoff
function calc_all_pair_forces_Tersoff!(calc::Calc)
    set_atom_forces_to_zero!(calc)
    for i in 1:calc.structure.noa
        calc.potential.precomputed = zeros(Float64,2, Int64(sum(@views calc.structure.neighbor_matrix[i,i+1:end])))
        calc_b_bp_full!(calc.potential,i)
        index = 1
        for ii in i+1:calc.structure.noa
            if calc.structure.neighbor_matrix[i,ii]
                partial_force = zeros(Float64, 3)
                b_ij = calc.potential.precomputed[1,index]
                bp_ij = calc.potential.precomputed[2,index]
                partial_force = calc_dU_i_dr_ij(calc.potential,i,ii,b_ij,bp_ij)   # dU_i_dr_ij
                calc.structure.atom_list[i].force += partial_force
                calc.structure.atom_list[ii].force -= partial_force
                index +=1
            end
        end
    end
end

## extra function for calculating the stress tensor
function calc_stress_tensor!(calc::Calc)
    calc.stress_tensor = zeros(Float64, 3,3)
    for ii in 2:calc.structure.noa
        if calc.structure.atom_list[ii].group == 0
            for i in 1:ii-1
                if calc.structure.neighbor_matrix[i,ii]
                    dist_vec = get_distance_between2atoms(calc.structure, i, ii)
                    calc.placeholder_force .= 0.0
                    calc_force!(calc, i, ii, stress=true)
                    calc.stress_tensor .+= dist_vec*calc.placeholder_force'
                end
            end
        end
    end
    calc.stress_tensor = (0.5)*calc.stress_tensor*(-1.0)/(calc.structure.box.lx*calc.structure.box.ly*calc.structure.box.lz)
end

#distribution function to the other calc force functions
# i  ... center atom number
# ii ... neighbor atom number
function calc_force!(calc::Calc, i::Int64, ii::Int64; stress::Bool=false)
    if calc.type_potential == 7 #Harm_ring_ring
        calc_force_hrr!(calc, i, ii)
    elseif calc.type_potential == 13 #Yukawa
        if stress || calc.structure.atom_list[i].norms[ii-i] < calc.potential.Rc
            calc_force_yuk!(calc, i, ii)
        end
    elseif calc.type_potential == 33 #LJ
        if stress || calc.structure.atom_list[i].norms[ii-i] < calc.potential.Rc
            calc_force_LJ!(calc, i, ii)
        end
    elseif calc.type_potential == 21 #Tersoff
        calc_force_ter!(calc.potential,i,ii)
    end
end

#this is the function for LJ since otherwise there are circular dependcys when giving the type of calc
function calc_force_LJ!(calc::Calc, i::Int64, ii::Int64)
    if calc.structure.atom_list[i].norms[ii-i] < calc.potential.Rc
        # current id of atom1 and atom2 -> find parameters
        pos = find_pos(calc.potential, calc.structure.atom_list[i].type, calc.structure.atom_list[ii].type)
        # calculate derivative
        dU_dr = calc_V_prime(calc.potential, calc.structure.atom_list[i].norms[ii-i], calc.potential.epsilon[pos], calc.potential.sigma[pos])
        # find shifting correction
        du_dr_corr = calc.potential.V_prime_Rc[pos]
        # calculate force
        calc.placeholder_force .= (du_dr_corr-dU_dr)/calc.structure.atom_list[i].norms[ii-i] .* @views calc.structure.atom_list[i].distances[:,ii-i]
    end
end

## calculate the inter-particle force
function calc_force_yuk!(calc::Calc, i::Int64, ii::Int64)
    dU_dr = 0.0
    if calc.structure.atom_list[i].norms[ii-i] < calc.potential.Rc
        # current id of atom1 and atom2 -> find parameters
        pos = find_pos(calc.potential, calc.structure.atom_list[i].type, calc.structure.atom_list[ii].type)
        # calculate derivative

        ## calculate the inverse of r
        dU_dr = calc_V_prime(calc.potential, calc.structure.atom_list[i].norms[ii-i], calc.potential.params[1][pos], calc.potential.params[2][pos])

        # find shifting correction
        du_dr_corr = calc.potential.V_prime_Rc[pos]
        # calculate force
        calc.placeholder_force .= (du_dr_corr- dU_dr)/calc.structure.atom_list[i].norms[ii-i] .* @views calc.structure.atom_list[i].distances[:, ii-i]
    end
end

## calculate the inter-particle force
function calc_force_hrr!(calc::Calc, i::Int64, ii::Int64)
    r_0 = calc.potential.R_0_mat[calc.structure.atom_list[i].type, calc.structure.atom_list[ii].type]
    dU_dr = -calc_V_prime(calc.potential, calc.structure.atom_list[i].norms[ii-i], r_0)
    @views calc.placeholder_force .= dU_dr/calc.structure.atom_list[i].norms[ii-i] .* @views calc.structure.atom_list[i].distances[:,ii-i]
end
