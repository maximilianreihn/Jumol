## Structure of the ensemble
######################################
######################################
######################################
#
#
# (c) Franz Bamer, Ivy Wu Mai-2022
######################################

include("atom.jl")
include("box.jl")
include("read.jl")
include("cell.jl")

mutable struct Structure
	## cutoff radius and skin
	rc::Float64
	rskin::Float64
	## periodic boundary conditions
	pbx::Int64
	pby::Int64
	pbz::Int64
	# linked cell boolean
	linked_cells_bool::Bool
	## load file reader
	reader::Read
	## number of atoms
	noa::Int64
	## atom list
	atom_list::Vector{Atom}
	## Box dimensions
	box::Box
	## constructor
	placeholder_pos::MVector{3,Float64} 
	placeholder_t::MVector{3,Float64}
	index_counter::Int64
	index_cache::Vector{Tuple{Int64, Int64}}
	total_mass::Float64
	Structure(;rc=10.0, rskin=1.0,
			  pbx=0, pby=0, pbz=0,
			  linked_cells_bool=true,
			  reader=Read(),
			  noa=0,
			  atom_list=Atom[],
			  box=Box(),
			  placeholder_pos=MVector{3,Float64}(zeros(Float64,3)),
			  placeholder_t=MVector{3,Float64}(zeros(Float64,3)),
			  index_counter=0,
			  index_cache=Tuple{Int64, Int64}[], total_mass=0.0) = new(rc, rskin,
											pbx, pby, pbz,
											linked_cells_bool,reader,noa,atom_list,box,
											placeholder_pos,placeholder_t,index_counter,index_cache,total_mass)
	neighbor_matrix::Matrix{Bool}
	types1::Int64
	cell_list::Vector{Cell}
	cell_neighbour_array::Matrix{Int64}
	reduced_coords::Vector{Float64}
	reduced_basis::Matrix{Float64}
	reduced_r_0::Vector{Float64}
	U_original::Float64
	reduced_dim::Int64
	all_grous_zero::Bool
end


##############################################
#### data loading and reading
##############################################
## load the box from a lammpstrj text file
function read_lammpstrj!(structure::Structure, filename::String; timestep_start=1, triclinic=false, include_vel=true)
	(structure)
	get_lines!(structure.reader, filename)
	read_lammpstrj_reader!(structure.reader, structure, structure.reader.line_start_list[timestep_start], triclinic=triclinic, include_vel=include_vel)
	set_box_basis_vectors!(structure.box)
end

## load the box from a lammps type box input file
function read_lammps_box!(structure::Structure,filename::String)
	(structure)
	read_box_lammps!(structure, filename)
end

## load an atom by hand
function add_atom_by_hand!(structure::Structure,
			number::Int64, type::Int64, x::Float64, 
			y::Float64, z::Float64, vx::Float64, vy::Float64, vz::Float64; group=0)
	push!(structure.atom_list, Atom(number,type))
	structure.atom_list[structure.noa+1].pos = [x,y,z]
	structure.atom_list[structure.noa+1].vel = [vx,vy,vz]
	structure.atom_list[structure.noa+1].acc = [0.0,0.0,0.0]
	structure.atom_list[structure.noa+1].group = group
	structure.noa += 1
end

## load the box dimensions by hand
function create_box_by_hand!(structure::Structure,
			lx::Float64, ly::Float64, lz::Float64,
			lxy::Float64, lyz::Float64, lxz::Float64)
	structure.box.lx = lx
	structure.box.ly = ly
	structure.box.lz = lz
	structure.box.lxy = lxy
	structure.box.lyz = lyz
	structure.box.lxz = lxz
	set_box_basis_vectors!(structure.box)
end

## get the positions in a matrix 3 x noa
function get_atom_pos_mat(structure::Structure)
	atom_pos_mat = zeros(Float64, structure.noa,3)
	for i in 1:structure.noa
		atom_pos_mat[i,:] .= @views structure.atom_list[i].pos[1:3]
	end
	return atom_pos_mat
end

## get global position vector noa*3 x 1
function get_global_atom_pos_vec!(structure::Structure, atom_pos_vec::Vector{Float64})
	for i in 1:structure.noa
		@views atom_pos_vec[(i-1)*3 + 1 : i*3] .= structure.atom_list[i].pos
	end
end
function get_global_atom_pos_vec(structure::Structure)
	atom_pos_vec = zeros(Float64, structure.noa*3)
	for i in 1:structure.noa
		@views atom_pos_vec[(i-1)*3 + 1 : i*3] .= structure.atom_list[i].pos
	end
	return atom_pos_vec
end

## set atomic positions according to a global position vector
function set_global_pos_vec!(structure::Structure, pos_vec::Vector{Float64})
	for i in 1:structure.noa
		@views structure.atom_list[i].pos .= pos_vec[(i-1)*3 + 1 : i*3]
	end
end

## get global force vector noa*3 x 1
function get_global_force_vec(structure::Structure, force_vec::Vector{Float64})
	for i in 1:structure.noa
		@views force_vec[(i-1)*3 + 1 : i*3] = structure.atom_list[i].force
	end
	return force_vec
end

## set atomic forces according to a global force vector
function set_global_force_vec!(structure::Structure, force_vec::Vector{Float64})
	for i in 1:structure.noa
		@views structure.atom_list[i].force .= force_vec[(i-1)*3 + 1 : i*3]
	end
end

##############################################
#### distance and neighbor list calculations
##############################################
## initialize structure (must be run once)
function initialize_structure!(structure::Structure)
	## create matrices for distances for every atom (4,noa-num_atom)
	#sort!(structure.atom_list, by=x->x.type) #first type 1 and then type 2
	#structure.types1 = 2*structure.noa -sum(Int64[a.type for a in structure.atom_list])
	"""
	structure.index_pairs = [filter(t -> t[1] != t[2] && t[1]<t[2], vec(collect(Iterators.product(1:structure.types1, 1:structure.types1)))), # all pairs of type 1  
							filter(t -> t[1] != t[2] && t[1]<t[2], vec(collect(Iterators.product(structure.types1+1:structure.noa, structure.types1+1:structure.noa)))), # all pairs of type 2
							vec(collect(Iterators.product(1:structure.types1, structure.types1+1:structure.noa)))] # all mixed pairs
	""";

	for i in 1:structure.noa
		structure.atom_list[i].distances = zeros(Float64, 3, structure.noa-i)
		structure.atom_list[i].norms = zeros(Float64, structure.noa-i)
	end
	if sum([a.group==0 for a in structure.atom_list]) == structure.noa
		structure.all_grous_zero = true
	else
		structure.all_grous_zero = false
	end
	structure.neighbor_matrix = zeros(Bool, structure.noa, structure.noa)
	## put atoms back to box if we have negative coordinates
	put_atoms_back_to_box!(structure)
	## initialize cell lists
	update_cell!(structure)
	construct_neigh_cell_list!(structure)
	linked_cell_list!(structure)
	update_distances!(structure)
end

## put atoms back to box if they leave the box
function put_atoms_back_to_box!(structure::Structure)
	for i in 1:structure.noa
		## to do (transformation into the triclinic coordinate system needed!!!!!!)
		if structure.atom_list[i].pos[1] > (structure.box.h1[1] +
			                                structure.box.h2[1]/structure.box.h2[2]* structure.atom_list[i].pos[2] +
											structure.box.h3[1]/structure.box.h3[3]*structure.atom_list[i].pos[3])
			structure.atom_list[i].pos -= structure.box.h1
		elseif structure.atom_list[i].pos[1] < (0.0 +
			                                    structure.box.h2[1]/structure.box.h2[2]* structure.atom_list[i].pos[2] +
			                                    structure.box.h3[1]/structure.box.h3[3]*structure.atom_list[i].pos[3])
			structure.atom_list[i].pos += structure.box.h1
		end
		if structure.atom_list[i].pos[2] > structure.box.h2[2]
			structure.atom_list[i].pos -= structure.box.h2
		elseif structure.atom_list[i].pos[2] < 0
			structure.atom_list[i].pos += structure.box.h2
		end
		if structure.atom_list[i].pos[3] > structure.box.h3[3]
			structure.atom_list[i].pos -= structure.box.h3
		elseif structure.atom_list[i].pos[3] < 0
			structure.atom_list[i].pos += structure.box.h3
		end
	end
end
