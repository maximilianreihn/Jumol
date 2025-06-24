## Distance calculations
######################################
######################################
######################################
# full distance updates
# verlet distance updates
# linked cell distance updates
#
# (c) Franz Bamer, Ivy Wu Mai-2022, Max Reihn Oct 2024
######################################

## calculate all distances for the whole model
function update_distances!(structure::Structure)
	structure.index_counter = 0
	structure.neighbor_matrix .= false
	for i in 1:structure.noa
		if structure.linked_cells_bool == true
			## calculate distance update using the linked cell lists
			for cell_id in @views structure.cell_neighbour_array[structure.atom_list[i].atom_cell_index_vec[4],:]
				if cell_id != 0
					for ii in structure.cell_list[cell_id].atoms_in_cell_list
						if i < ii
							calc_distance_betw2atoms!(structure, i, ii)
							structure.neighbor_matrix[i, ii] = structure.neighbor_matrix[ii, i] = structure.atom_list[i].norms[ii-i] < structure.rc+structure.rskin
							if structure.neighbor_matrix[i, ii]
								structure.index_counter += 1
								if length(structure.index_cache) < structure.index_counter
									push!(structure.index_cache, (i,ii))
								else
									structure.index_cache[structure.index_counter] = (i,ii)
								end
							end
						end
					end
				end
			end
		else
			## calculate full distance update
			for ii in i+1:structure.noa
				calc_distance_betw2atoms!(structure, i, ii)
				structure.neighbor_matrix[i, ii] = structure.neighbor_matrix[ii, i] = structure.atom_list[i].norms[ii-i] < structure.rc+structure.rskin
				if structure.neighbor_matrix[i, ii]
					structure.index_counter += 1
					if length(structure.index_cache) < structure.index_counter
						push!(structure.index_cache, (i,ii))
					else
						structure.index_cache[structure.index_counter] = (i,ii)
					end
				end
			end
		end
	end
end

function update_index_cache(structure::Structure)
	structure.index_counter = 0
	for i in 1:structure.noa
		for ii in i+1:structure.noa
			if structure.neighbor_matrix[i, ii]
				structure.index_counter += 1
				if length(structure.index_cache) < structure.index_counter
					push!(structure.index_cache, (i,ii))
				else
					structure.index_cache[structure.index_counter] = (i,ii)
				end
			end
		end
	end
end

## update distances in within the cutoff radius
function update_distances_partly!(structure::Jumol.Structure)
	for (i,ii) in @views structure.index_cache[1:structure.index_counter]
		calc_distance_betw2atoms!(structure, i, ii)
	end
end

## update the cell list and information
function update_cell!(structure::Structure)
	# modify the cutoff radius for shear condition
	cos_gamma = sqrt(structure.box.lxy * structure.box.lxy + structure.box.ly * structure.box.ly) / structure.box.ly
	r_cutoff = (structure.rc+structure.rskin)/cos_gamma
	# subdivide global box into n smaller cells
	n_cells = max.(Vector{Int64}([floor(Int64, structure.box.lx/r_cutoff),
		floor(Int64, structure.box.ly/r_cutoff),
		floor(Int64, structure.box.lz/r_cutoff)]),
		ones(Int64, 3))
	# erase the original cell list
	structure.cell_list = Vector{Cell}([Cell(cell_id) for cell_id in 1:prod(n_cells)])
end

#for tmap 
function is_in_box!(structure::Structure, i::Int64, upper_limits::Vector{Int64})
	for j in 1:3
		if structure.atom_list[i].atom_cell_index_vec[j] < 0
			structure.atom_list[i].atom_cell_index_vec[j] = upper_limits[j]-1
		elseif structure.atom_list[i].atom_cell_index_vec[j] >= upper_limits[j]
			structure.atom_list[i].atom_cell_index_vec[j] = 0
		end
	end
end

## create list of linked cells
function linked_cell_list!(structure::Structure)
	# modify the cutoff radius for shear condition
	tan_gamma = structure.box.lxy/structure.box.ly
	cos_gamma = sqrt(structure.box.lxy^2 + structure.box.ly^2) / structure.box.ly
	r_cutoff = (structure.rc+structure.rskin)/cos_gamma
	# subdivide global box into n smaller cells
	n_cells = max.(Vector{Int64}([floor(Int64, structure.box.lx/r_cutoff),
			floor(Int64, structure.box.ly/r_cutoff),
			floor(Int64, structure.box.lz/r_cutoff)]),
			ones(Int64, 3))
	# cell length  !Fixme! if 2d, there will be a zero in n_cell_z
	l_cells = Float64[geo/n_cells[i] for (i ,geo) in enumerate([structure.box.lx,structure.box.ly,structure.box.lz])]
	#Threads.@threads actually this slows down the process also the pushing needs to be changed if done with threads
	map(a -> a.pos[1]-=tan_gamma * a.pos[2], structure.atom_list) 
	map(a -> a.atom_cell_index_vec = Vector{Int64}(undef, 4), structure.atom_list)
	#all_max = maximum([norm(a.pos) for a in structure.atom_list])
	#println("Here we have ", l_cells, "all max ", all_max)
	map(a -> a.atom_cell_index_vec[1:3] .= floor.(Int64, a.pos./l_cells), structure.atom_list)
	for i in 1:structure.noa
		is_in_box!(structure, i, n_cells)
		cell_index = structure.atom_list[i].atom_cell_index_vec[1]*n_cells[2]*n_cells[3] + 
					structure.atom_list[i].atom_cell_index_vec[2]*n_cells[3] + 
					structure.atom_list[i].atom_cell_index_vec[3] + 1
		#println("Cell index is ", cell_index,' ', structure.atom_list[i].pos, ' ', l_cells)
		structure.atom_list[i].atom_cell_index_vec[4] = Int64(cell_index)
		push!(structure.cell_list[Int64(cell_index)].atoms_in_cell_list, i)
	end
	map(a -> a.pos[1]+=tan_gamma * a.pos[2], structure.atom_list)
end

## create list of neighboring cells
function construct_neigh_cell_list!(structure::Structure)
	# modify the cutoff radius for shear condition
	cos_gamma = sqrt(structure.box.lxy * structure.box.lxy + structure.box.ly * structure.box.ly) / structure.box.ly
	r_cutoff = (structure.rc+structure.rskin)/cos_gamma
	# subdivide global box into n smaller cells
	n_cells = Vector{Int64}([floor(Int64, structure.box.lx/r_cutoff), floor(Int64, structure.box.ly/r_cutoff), floor(Int64, structure.box.lz/r_cutoff)])
	n_cells .= map(x->x==0 ? x=1 : x=x, n_cells)
	l_cells = Float64[geo/n_cells[i] for (i ,geo) in enumerate([structure.box.lx,structure.box.ly,structure.box.lz])]
	n_cell_xyz = prod(n_cells)
	#
    # construct lists of neighbouring cells
	#println("Here we habe n cells ", n_cells, " ",structure.box.lx," ", structure.box.ly," ",structure.box.lz," ",r_cutoff)
	neighboring_cell_array = zeros(Int64, Int64(n_cell_xyz), 27)
	structure.cell_neighbour_array = zeros(Int64, Int64(n_cell_xyz), 27)
	# current cell
	for i in 0:n_cells[1]-1
		for ii in 0:n_cells[2]-1
			for iii in 0:n_cells[3]-1
				c = Int64(i*n_cells[2]*n_cells[3] + ii*n_cells[3] + iii + 1)
				# neighbour cells
				count_neigh = 1
				for j in i-1:i+1
					for jj in ii-1:ii+1
						for jjj in iii-1:iii+1
							#
							if structure.pbx == 1
								if j<0
									j += n_cells[1]
								elseif j==n_cells[1]
									j = 0
								end
							end
							#
							if structure.pby == 1
								if jj<0
									jj += n_cells[2]
								elseif jj==n_cells[2]
									jj = 0
								end
							end
							#
							if structure.pbz == 1
								if jjj<0
									jjj += n_cells[3]
								elseif jjj==n_cells[3]
									jjj = 0
								end
							end
							#
							if j>=0 && jj>=0 && jjj>=0 && j<n_cells[1] && jj<n_cells[2] && jjj<n_cells[3]
								c1 = Int64(j*n_cells[2]*n_cells[3] + jj*n_cells[3] + jjj +1)
								neighboring_cell_array[c,count_neigh] = c1
								structure.cell_neighbour_array[c,count_neigh] = c1
							end
							count_neigh += 1
							#
						end
					end
				end
			end
		end
	end
end


## calculate the distance between two atoms
function calc_distance_betw2atoms!(structure::Structure, num1::Int64, num2::Int64)
	structure.placeholder_pos .= structure.atom_list[num2].pos
	dot_prods = 0.0
	structure.placeholder_t .= structure.atom_list[num1].pos .- structure.placeholder_pos
	## periodic boundary conditions in e1
	if structure.pbx == 1
		dot_prods = -2.0*dot(structure.placeholder_t, structure.box.e1)
		if dot_prods > structure.box.l1
			structure.placeholder_pos -= structure.box.h1
			structure.placeholder_t += structure.box.h1
		elseif dot_prods < -structure.box.l1
			structure.placeholder_pos += structure.box.h1
			structure.placeholder_t -= structure.box.h1
		end
	end
	## periodic boundary conditions in e2
	if structure.pby == 1
		dot_prods = -2.0*dot(structure.placeholder_t, structure.box.e2)
		if dot_prods > structure.box.l2
			structure.placeholder_pos -= structure.box.h2
			structure.placeholder_t += structure.box.h2
		elseif dot_prods < -structure.box.l2
			structure.placeholder_pos += structure.box.h2
			structure.placeholder_t -= structure.box.h2
		end
	end
	## periodic boundary conditions in e3
	if structure.pbz == 1
		dot_prods = -2.0*dot(t_2_1, structure.box.e3)
		if dot_prods > structure.box.l3
			structure.placeholder_pos -= structure.box.h3
			structure.placeholder_t += structure.box.h3
		end
		if dot_prods < -structure.box.l3
			structure.placeholder_pos += structure.box.h3
			structure.placeholder_t -= structure.box.h3
		end
	end
	@views structure.atom_list[num1].distances[:, num2-num1] .= structure.placeholder_t
	@views structure.atom_list[num1].norms[num2-num1] = norm(structure.placeholder_t)
end

## get the distance between two atoms
function get_distance_between2atoms(structure::Structure, num1::Int64, num2::Int64; norm_only::Bool = false)
	if norm_only
		return structure.atom_list[num1].norms[num2-num1]
	elseif num1 <= num2
		return @views structure.atom_list[num1].distances[:,num2-num1]
	else
		return @views -structure.atom_list[num2].distances[:,num1-num2]
	end
end
