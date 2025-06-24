## time integrators
######################################
######################################
######################################

using Random, Distributions

mutable struct Integrator
    calc::Calc
	units::Units
    delta_t::Float64
	output_filename::String
	steps_save::Int64
    Integrator(calc, units; 
			delta_t=1.0e-2,
		  	output_filename="None",
		  	steps_save=1) = new(calc,units,
										delta_t,
										output_filename,
										steps_save)
	## time step vector
	timestep_vec::Vector{Float64}
	## NVT vector
	temp_vec_target::Vector{Float64}
	temp_vec::Vector{Float64}
	zeta_vec::Vector{Float64}
end


## NVE ensemble
function run_nve!(int::Integrator, num_steps::Int64; quiet = false)
	## write file if necessary
	if int.output_filename != "None"
		open_file!(int.calc.structure.reader, int.output_filename, "w")
	end
	## starting NVE integration
	if !quiet
		println("####")
		println("starting NVE integration...")
	end
	## calculate initial force
	calc_all_pair_forces!(int.calc)
	## get storage for the calculation parameters
	# velocity at the half step
	vel_half = zeros(Float64, int.calc.structure.noa,3)
	# maximum walking distance
	max_dist = 0
	# maximum vel at each time step
	max_vel = 0
 	for i in 1:num_steps
		if !quiet
			print("\e[2K")
			print("\rstep: ", i)
		end
		# update position and calculate intermediate velocity
        for ii in 1:int.calc.structure.noa
            vel_half[ii,:] = int.calc.structure.atom_list[ii].vel + 0.5 * int.calc.structure.atom_list[ii].force/int.calc.structure.atom_list[ii].mass * int.units.force_factor*int.delta_t
            int.calc.structure.atom_list[ii].pos .+= vel_half[ii,:] * int.delta_t
        end
        # update distances
		for ii in 1:int.calc.structure.noa
			v = norm(int.calc.structure.atom_list[ii].vel)
			if v > max_vel
				max_vel = v
			end
		end
		max_dist += max_vel * int.delta_t
		max_vel = 0

		if max_dist > int.calc.structure.rskin * 0.5
			if int.calc.structure.linked_cells_bool == true
				update_cell!(int.calc.structure)
				linked_cell_list!(int.calc.structure)
				construct_neigh_cell_list!(int.calc.structure)
			end
	    	update_distances!(int.calc.structure)
			max_dist = 0
		else
	    	update_distances_partly!(int.calc.structure)
		end
        # update force
        calc_all_pair_forces!(int.calc)
        # update velocity
        for ii in 1:int.calc.structure.noa
	    	int.calc.structure.atom_list[ii].vel = vel_half[ii,:] + 0.5 * int.calc.structure.atom_list[ii].force/int.calc.structure.atom_list[ii].mass * int.units.force_factor* int.delta_t
		end
		put_atoms_back_to_box!(int.calc.structure)
		if int.output_filename != "None"
			if i%int.steps_save == 0
				Jumol.write_box!(int.calc.structure.reader,int.calc.structure,i)
			end
		end
    end
	## close output file
	if int.output_filename != "None"
		close_file!(int.calc.structure.reader)
	end
	return vel_half
end

## NVT ensemble
function run_nvt!(int::Integrator, num_steps::Int64, t_beg::Float64, t_end::Float64, tau::Float64; order_of_problem=2, zeta=0.0, quiet = false)
	## write file if necessary
	if int.output_filename != "None"
		open_file!(int.calc.structure.reader, int.output_filename, "w")
	end
	## set atom velocities
	set_atom_velocities(int, 2, t_beg)
	## get temperature history
	int.temp_vec_target= get_temp_vec(int,t_beg,t_end,num_steps)
	int.timestep_vec = LinRange(1,num_steps,num_steps)
	## starting NVT integration
	if !quiet
		println("####")
		println("starting NVT integration...")
	end
	## calculate initial force
	calc_all_pair_forces!(int.calc)
	## get storage for the calculation parameters
	# velocity at the half step
	vel_half = zeros(Float64, int.calc.structure.noa,3)
	# maximum walking distance
	max_dist = 0
	# maximum vel at each time step
	max_vel = 0
	# target temperature (must be defined as input parameter)
	T_tar = 0.0
	# initialize actual temperature
	T = 0.0
	# initialize kinetic energy
	KE = 0
	# degrees of freedom for 2d
	dof = order_of_problem * int.calc.structure.noa
	# temperature vector at all time steps
	int.temp_vec = zeros(Float64, num_steps)
	int.zeta_vec = zeros(Float64, num_steps)
	for i in 1:num_steps
		## temperature
		T_tar = int.temp_vec_target[i]
		if !quiet
			print("\e[2K")
			print("\rstep: ", i)
		end
		#update position and calculate intermediate velocity
    	for ii in 1:int.calc.structure.noa
        	vel_half[ii,:] = int.calc.structure.atom_list[ii].vel + 0.5 * (int.calc.structure.atom_list[ii].force /int.calc.structure.atom_list[ii].mass* int.units.force_factor - zeta * int.calc.structure.atom_list[ii].vel)* int.delta_t
			int.calc.structure.atom_list[ii].pos .+= vel_half[ii,:] * int.delta_t
		end
		# update distances
		for ii in 1:int.calc.structure.noa
			v = norm(int.calc.structure.atom_list[ii].vel)
			if v > max_vel
				max_vel = v
			end
		end
		max_dist += max_vel * int.delta_t
		max_vel = 0
		if max_dist > int.calc.structure.rskin*0.5
			if int.calc.structure.linked_cells_bool == true
                update_cell!(int.calc.structure)
            	linked_cell_list!(int.calc.structure)
            	construct_neigh_cell_list!(int.calc.structure)
            end
    		update_distances!(int.calc.structure)
			max_dist = 0
		else
    		update_distances_partly!(int.calc.structure)
		end
    	#update force
    	calc_all_pair_forces!(int.calc)
		# calculate kinetic energy
		KE = 0
		for ii in 1:int.calc.structure.noa
			KE += 0.5 * int.calc.structure.atom_list[ii].mass * dot(int.calc.structure.atom_list[ii].vel,int.calc.structure.atom_list[ii].vel)
		end
		# calculate temperature
		T = 2 * KE / (dof * int.units.Kb )
		zeta += 0.5 * int.delta_t * (T/T_tar - 1.0)/tau^2
		KE = 0
		for ii in 1:int.calc.structure.noa
			KE += 0.5 * int.calc.structure.atom_list[ii].mass  * dot(vel_half[ii,:],vel_half[ii,:])
		end
		T = 2*KE / (dof * int.units.Kb)
		zeta += 0.5 * int.delta_t * (T/T_tar - 1.0)/tau^2
    	#Update velocity
    	for ii in 1:int.calc.structure.noa
    		int.calc.structure.atom_list[ii].vel = (vel_half[ii,:] + 0.5 * int.calc.structure.atom_list[ii].force/int.calc.structure.atom_list[ii].mass* int.units.force_factor * int.delta_t )/(1 + 0.5 * int.delta_t * zeta )
		end
		int.temp_vec[i] = T
		int.zeta_vec[i] = zeta
		put_atoms_back_to_box!(int.calc.structure)
		if int.output_filename != "None"
			if i%int.steps_save == 0
				Jumol.write_box!(int.calc.structure.reader,int.calc.structure,i)
			end
		end
	end
	## close output file
	if int.output_filename != "None"
		close_file!(int.calc.structure.reader)
	end
	print("\n")
end


function run_npt(int::Integrator, num_steps::Int64)
	println("... npt not yet implemented ...")
end


function get_temp_vec(int::Integrator, t_beg::Float64, t_end::Float64, num_steps::Int64)
	temp_target_vec = zeros(Float64, num_steps)
	for i in 1:num_steps
		temp_target_vec[i] = t_beg + (t_end-t_beg)/Float64(num_steps) * (i-1)
	end
	return temp_target_vec
end

## Bolzmann distribution
function set_atom_velocities(int::Integrator, order_of_problem::Int64, T::Float64)
	for atom in int.calc.structure.atom_list
		mu_vel = 0.0
		sigma_vel = sqrt(int.units.Kb * T/atom.mass)
		dist = Normal(mu_vel,sigma_vel)
		if order_of_problem == 2
			atom.vel = [rand(dist), rand(dist), 0.0]
		elseif order_of_problem == 3
			atom.vel = [rand(dist), rand(dist), rand(dist)]
		end
	end
end
