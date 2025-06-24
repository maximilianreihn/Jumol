## calculation of Hessian and related quantities

mutable struct Hessian
    calc::Calc
    Hessian(calc) = new(calc)
    hessian_mat::Matrix{Float64}
    evals::Vector{Float64}
    evecs::Matrix{Float64}
    hessian_mat_sparse::Matrix{Float64}
end

function calc_hessian_anal!(hessian::Hessian)
    hessian.hessian_mat = zeros(Float64, 2*hessian.calc.structure.noa, 2*hessian.calc.structure.noa)
    temp_mat = @MMatrix zeros(Float64, 2, 2)  # StaticArray for fixed-size matrices
    temp_diag_mat = @MMatrix zeros(Float64, 2, 2)
    eyes =SMatrix{2,2}(Matrix{Float64}(I,2,2))
    for i in 1:hessian.calc.structure.noa
        temp_diag_mat .= 0.0
        for j in i+1:hessian.calc.structure.noa
            if hessian.calc.structure.neighbor_matrix[i,j]
                norm_val = hessian.calc.structure.atom_list[i].norms[j-i]
                norm_fac = 1 / norm_val
                pos = find_pos(hessian.calc.potential, hessian.calc.structure.atom_list[i].type, hessian.calc.structure.atom_list[j].type)
                eps = hessian.calc.potential.epsilon[pos]
                sig = hessian.calc.potential.sigma[pos]

                V_prime = calc_V_prime(hessian.calc.potential, norm_val, eps, sig)
                V_pp = calc_V_pp(norm_val, eps, sig)
                factor = -norm_fac * V_prime
                term2 = -(norm_fac^2 * V_pp - norm_fac^3 * V_prime)
                dist_vec = @view hessian.calc.structure.atom_list[i].distances[1:2, j-i]
                temp_mat .= factor * eyes + term2 * (dist_vec * dist_vec')
                @views hessian.hessian_mat[2*i-1 : 2*i, 2*j-1 : 2*j] .= temp_mat  # Upper diagonal
                @views hessian.hessian_mat[2*j-1 : 2*j, 2*i-1 : 2*i] .= temp_mat'  # Symmetric entries
                temp_diag_mat .-= temp_mat  # Cache diagonal entries
                @views hessian.hessian_mat[2*j-1 : 2*j, 2*j-1 : 2*j] .-= temp_mat  # Diagonal adjustment
            end
        end
        @views hessian.hessian_mat[2*i-1 : 2*i, 2*i-1 : 2*i] .+= temp_diag_mat  # Update diagonal entries
    end
end

function calc_hessian_anal_3D!(hessian::Hessian)
    hessian.hessian_mat = zeros(Float64, 3*hessian.calc.structure.noa, 3*hessian.calc.structure.noa)
    temp_mat = @MMatrix zeros(Float64, 3, 3)  # StaticArray for fixed-size matrices
    temp_diag_mat = @MMatrix zeros(Float64, 3, 3)
    eyes =SMatrix{3,3}(Matrix{Float64}(I,3,3))
    for i in 1:hessian.calc.structure.noa
        temp_diag_mat .= 0.0
        for j in i+1:hessian.calc.structure.noa
            if hessian.calc.structure.neighbor_matrix[i,j]
                norm_val = hessian.calc.structure.atom_list[i].norms[j-i]
                norm_fac = 1 / norm_val
                pos = find_pos(hessian.calc.potential, hessian.calc.structure.atom_list[i].type, hessian.calc.structure.atom_list[j].type)
                eps = hessian.calc.potential.epsilon[pos]
                sig = hessian.calc.potential.sigma[pos]

                V_prime = calc_V_prime(hessian.calc.potential, norm_val, eps, sig)
                V_pp = calc_V_pp(norm_val, eps, sig)
                factor = -norm_fac * V_prime
                term2 = -(norm_fac^2 * V_pp - norm_fac^3 * V_prime)
                dist_vec = @view hessian.calc.structure.atom_list[i].distances[1:3, j-i]
                temp_mat .= factor * eyes + term2 * (dist_vec * dist_vec')
                @views hessian.hessian_mat[3*i-1 : 3*i, 3*j-1 : 3*j] .= temp_mat  # Upper diagonal
                @views hessian.hessian_mat[3*j-1 : 3*j, 3*i-1 : 3*i] .= temp_mat'  # Symmetric entries
                temp_diag_mat .-= temp_mat  # Cache diagonal entries
                @views hessian.hessian_mat[3*j-1 : 3*j, 3*j-1 : 3*j] .-= temp_mat  # Diagonal adjustment
            end
        end
        @views hessian.hessian_mat[3*i-1 : 3*i, 3*i-1 : 3*i] .+= temp_diag_mat  # Update diagonal entries
    end
end

## main function to calculate Hessian
function calc_hessian!(hessian::Hessian; delta::Float64=1.0e-5)

    hessian.hessian_mat = zeros(Float64, 2*hessian.calc.structure.noa,2*hessian.calc.structure.noa)

    pos2 = @MVector zeros(Float64, 3)
    t_2_1 = @MVector zeros(Float64, 3)

    for col in 1:hessian.calc.structure.noa  # row -> pos_atom
        for row in 1:col                # col -> force_atom
            if row == col || hessian.calc.structure.neighbor_matrix[row, col]
                calc_hessian_block_matrix!(hessian, pos2, t_2_1, row, col, delta)
            end
        end
    end
    ## symmetric entries
    for col in 1:hessian.calc.structure.noa*2
        for row in (col+1):hessian.calc.structure.noa*2
            hessian.hessian_mat[row,col] = hessian.hessian_mat[col,row]
        end
    end
    hessian.hessian_mat .*= -1.0
end

## calculate hessian block matrices
function calc_hessian_block_matrix!(hessian::Hessian, pos2::MVector{3, Float64}, t_2_1::MVector{3, Float64}, row::Int64, col::Int64, delta::Float64)

    # F_x, pos direction-x (11),  F_y, pos direction-x(21)
    calc_force_derivative_hessian!(hessian, pos2, t_2_1, row, 1, col, delta)
    # F_x, pos direction-y(12) ,  F_y, pos direction-y(22)
    calc_force_derivative_hessian!(hessian, pos2, t_2_1, row, 2, col, delta)
end

## function for calculating force derivative
function calc_force_derivative_hessian!(hessian::Hessian, pos2::MVector{3, Float64}, t_2_1::MVector{3, Float64}, pos_atom_index::Int64, pos_dir::Int64, force_atom_index::Int64, delta::Float64)
    ## move atom in positive direction by delta
    hessian.calc.structure.atom_list[pos_atom_index].pos[pos_dir] += delta
    update_distances_hessian!(hessian, pos2, t_2_1, force_atom_index, pos_atom_index)
    hessian.calc.structure.atom_list[force_atom_index].force .= 0.0 #zeros(Float64, 3) # set_atom_forces_to_zero for "force_atom"
    calc_atom_force_hessian!(hessian,force_atom_index, pos_atom_index)
    f_dir_plus = hessian.calc.structure.atom_list[force_atom_index].force #[force_dir]

    ## move atom in negative direction by delta
    hessian.calc.structure.atom_list[pos_atom_index].pos[pos_dir] -= 2.0 * delta
    update_distances_hessian!(hessian, pos2, t_2_1,force_atom_index, pos_atom_index)
    hessian.calc.structure.atom_list[force_atom_index].force = @MVector zeros(Float64, 3) # set_atom_forces_to_zero for "force_atom"
    calc_atom_force_hessian!(hessian,force_atom_index, pos_atom_index)
    f_dir_minus = hessian.calc.structure.atom_list[force_atom_index].force #[force_dir]

    ## move atom back to its initial position and re-calcuate initial force (why?)
    hessian.calc.structure.atom_list[pos_atom_index].pos[pos_dir] += delta
    ## may_be_update distances
    update_distances_hessian!(hessian, pos2, t_2_1,force_atom_index, pos_atom_index)

    f_derivative = (f_dir_plus[1:2]-f_dir_minus[1:2])./(2.0 * delta)

    hessian.hessian_mat[2*pos_atom_index-1:2*pos_atom_index,2*force_atom_index-2+pos_dir] = f_derivative
end

## function for calculating force on atom with index i
function calc_atom_force_hessian!(hessian::Hessian, i::Int64, j::Int64)
    if hessian.calc.type_potential == 21
        cal_three_body_atom_force_hessian!(hessian, i)
    else
        cal_two_body_atom_force_hessian!(hessian, i)
        #cal_two_body_atom_force_hessian_optimized!(hessian::Hessian, i, j)
    end
end

## calculate two body atom force optimized
function cal_two_body_atom_force_hessian_optimized!(hessian::Hessian, i::Int64, j::Int64)
    force = zeros(Float64, 3)
    if j == i
        for ii in i+1:hessian.calc.structure.noa
            if hessian.calc.structure.neighbor_matrix[i,ii]
                force .= 0.0
                calc_force!(force, hessian.calc, i, ii)
                hessian.calc.structure.atom_list[i].force .+= force
            end
        end
    else
        force .= 0.0
        calc_force!(force, hessian.calc, min(i,j), max(i,j))
        hessian.calc.structure.atom_list[i].force .+= force
    end
end
## function for calculating two-body force on atom with index 'i'
function cal_two_body_atom_force_hessian!(hessian::Hessian, i::Int64)
    force = zeros(Float64, 3)
    for ii in 1:hessian.calc.structure.noa
        if hessian.calc.structure.neighbor_matrix[i,ii]
            force .= 0.0
            if i < ii
                calc_force!(force, hessian.calc,i,ii)
                hessian.calc.structure.atom_list[i].force .+= force
            else
                calc_force!(force, hessian.calc,ii,i)
                hessian.calc.structure.atom_list[i].force .-= force
            end
        end
    end
end


## function for calculating three-body force (Tersoff) on atom with index i
function cal_three_body_atom_force_hessian!(hessian::Hessian, i::Int64)
    if hessian.calc.type_potential == 21
        ## calculate the sum of terms -> Σ dU_i_dr_ij
        size_ref = sum(hessian.calc.structure.neighbor_matrix[i, i+1:end])
        hessian.calc.potential.precomputed = zeros(Float64, 2,size_ref)
        calc_b_bp_full!(hessian.calc.potential,i)
        counter = 1
        for ii in i+1:hessian.calc.structure.noa
            if hessian.calc.structure.neighbor_matrix[i,ii]  # ii = j
                partial_force = zeros(Float64, 3)
                b_ij = hessian.calc.potential.precomputed[1,counter]
                bp_ij = hessian.calc.potential.precomputed[2,counter]
                partial_force = calc_dU_i_dr_ij(hessian.calc.potential,i,ii,b_ij,bp_ij)   # dU_i_dr_ij term
                hessian.calc.structure.atom_list[i].force .+= partial_force
                counter +=1
            end
        end
        ## calculate the sum of terms -> -Σ dU_j_dr_ji

        for ii in i+1:hessian.calc.structure.noa
            if hessian.calc.structure.neighbor_matrix[i,ii]
                size_ref = sum(hessian.calc.structure.neighbor_matrix[ii, ii+1:end])
                hessian.calc.potential.precomputed = zeros(Float64, 2, size_ref)
                calc_b_bp_full!(hessian.calc.potential,ii)
                counter = sum(hessian.calc.structure.neighbor_matrix[i, i+1:ii])
                partial_force = zeros(Float64, 3)
                b_ij = hessian.calc.potential.precomputed[1,counter]
                bp_ij = hessian.calc.potential.precomputed[2,counter]
                partial_force = calc_dU_i_dr_ij(hessian.calc.potential,ii,i,b_ij,bp_ij)   # dU_j_dr_ji term
                hessian.calc.structure.atom_list[i].force .-= partial_force
            end
        end
    end
end

#### eigenvalue eigenvector calcuation
## function for calculating eigen values and eigen vectors
function solve_eig!(hessian::Hessian; calc_eigen_vector = false, num_eigen_vectors = -1)
    if num_eigen_vectors == -1
        ef = eigen(Symmetric(hessian.hessian_mat))  # all eigen-values/eigen-vectors
    elseif num_eigen_vectors > 0
        ef = eigen(Symmetric(hessian.hessian_mat), 1:num_eigen_vectors)   #k smallest eigenvalues/vectors
    else
        print("Error while calling solve_eig")
    end
    hessian.evals = ef.values
    if calc_eigen_vector == true
        hessian.evecs = ef.vectors
    end
end

## function for calculating eigen values and eigen vectors of the dynamical matrix
function solve_eig_dynamical_mat(hessian::Hessian, inv_mass_mat::Array{Float64}; calc_eigen_vector = false, num_eigen_vectors = -1)
    D = inv_mass_mat * hessian.hessian_mat
    evals_D = nothing
    evecs_D = nothing
    if num_eigen_vectors == -1
        ef = eigen(Symmetric(D))  # all eigen-values/eigen-vectors
    elseif num_eigen_vectors > 0
        ef = eigen(Symmetric(D), 1:num_eigen_vectors)   #k smallest eigenvalues/vectors
    else
        print("Error while calling solve_eig of the dynamical matrix")
    end
    evals_D = ef.values
    if calc_eigen_vector == true
        evecs_D = ef.vectors
    end
    return evals_D, evecs_D
end

## update distance functions specific to two body force Hessian calculation
function update_distances_hessian!(hessian::Hessian, pos2::MVector{3, Float64}, t_2_1::MVector{3, Float64}, force_atom_index::Int64, pos_atom_index::Int64)
    if hessian.calc.type_potential == 21
        # update_distances(hessian.calc.structure)
        update_distances_partly!(hessian.calc.structure)
    else
        # update_distances(hessian.calc.structure)
        updates_distances_hessian_two_body_force!(hessian.calc.structure, pos2, t_2_1, force_atom_index, pos_atom_index)
    end
end

function updates_distances_hessian_two_body_force!(structure::Structure, pos2::MVector{3, Float64}, t_2_1::MVector{3, Float64}, i::Int64, j::Int64)
    if j == i
        for ii in 1:structure.noa
            if structure.neighbor_matrix[i,ii]
                update_distances_between_two_atoms!(structure, pos2, t_2_1, min(i,ii), max(i,ii))
            end
        end
    else
        update_distances_between_two_atoms!(structure, pos2, t_2_1, i, j)
    end
end

function update_distances_between_two_atoms!(structure::Structure, pos2::MVector{3, Float64}, t_2_1::MVector{3, Float64}, i::Int64, ii::Int64)
	if ii > i
		calc_distance_betw2atoms!(structure, pos2, t_2_1, i, ii)
    else
        calc_distance_betw2atoms!(structure, pos2, t_2_1, ii, i)
	end
end
