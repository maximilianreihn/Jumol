## Atom class
######################################
######################################
######################################
#import Core.Array


mutable struct Atom
	number::Int64
	type::Int64
	group::Int64
	calc_stress::Bool
	mass::Float64
	pot::Float64
	force::MVector{3,Float64}
	Atom(number,type;
		group=0,
		calc_stress=true,
		mass = 1.0,
		pot = 0.0,
		force = MVector{3,Float64}(zeros(Float64,3))) = new(number,type,group,calc_stress,mass,pot,force)

	pos::MVector{3,Float64}
	vel::MVector{3,Float64}
	acc::MVector{3,Float64}
	##
	distances::Matrix{Float64}
	norms::Vector{Float64}
	##
	## cell index, to which atom belongs to
	atom_cell_index_vec::Vector{Int64}
end