## Box class
######################################
######################################
######################################


mutable struct Box
	lx::Float64
	ly::Float64
	lz::Float64
	lxy::Float64
	lyz::Float64
	lxz::Float64
	l1::Float64
	l2::Float64
	l3::Float64
	## box dimensions
	## basis vectors
	h1::MVector{3, Float64}
	h2::MVector{3, Float64}
	h3::MVector{3, Float64}
	## length of the tricilic box
	## unit vectors
	e1::MVector{3, Float64}
	e2::MVector{3, Float64}
	e3::MVector{3, Float64}
	function Box(;lx=0.,ly=0.,lz=0.,lxy=0.,lyz=0.,lxz=0.,l1=0.,l2=0.,l3=0.) 
		h1 = MVector{3,Float64}([lx, 0., 0.])
		h2 = MVector{3,Float64}([lxy, ly, 0.])
		h3 = MVector{3,Float64}([lxz, lyz, lz])
		e1 = MVector{3,Float64}(zeros(Float64, 3))
		e2 = MVector{3,Float64}(zeros(Float64, 3))
		e3 = MVector{3,Float64}(zeros(Float64, 3))
		return new(lx,ly,lz,lxy,lyz,lxz,l1,l2,l3,h1,h2,h3,e1,e2,e3)
	end
end

function set_box_basis_vectors!(box::Box)
	## box basis vectors
	box.h1[1] = box.lx
	box.h2[1] = box.lxy
	box.h2[2] = box.ly
	box.h3[1] = box.lxz
	box.h3[2] = box.lyz
	box.h3[3] = box.lz
	## lenght of the triclinic box
	box.l1 = norm(box.h1)
	box.l2 = norm(box.h2)
	box.l3 = norm(box.h3)
	## box unit vectors
	box.e1 .= box.h1 ./ box.l1
	box.e2 .= box.h2 ./ box.l2
	box.e3 .= box.h3 ./ box.l3
end
