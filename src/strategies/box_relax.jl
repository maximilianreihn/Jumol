## Box relax
######################################
######################################
######################################
# relaxation of the box
#    - volumetric
#    - biaxial
#    - simple_shear
######################################

mutable struct Box_relax
    molcalc::Calc
    deform::Aff_deform
    max_num_iter::Int64
    energy_tol::Float64
    stress_tol::Float64
    Box_relax(molcalc,deform; max_num_iter=100,
              energy_tol=1.0e-5,stress_tol=1.0e-8) = new(molcalc,deform,
                                                         max_num_iter,energy_tol,
                                                         stress_tol)
end

## calculate pressure
function get_press(box_relax::Box_relax)
    press = 0.0
    for i in 1:3
        press += box_relax.molcalc.stress_tensor[i,i]
    end
    return press/3.0
end

function relax_vol!(box_relax::Box_relax, delta::Float64)
    ## initial direction of relaxation
    press = get_press(box_relax)
    direc = sign(press*(-1.0))
    for ii in 1:box_relax.max_num_iter
        println("######################")
        println("volumentric relaxation")
        ## get stress tensor
        Jumol.calc_all_pair_forces!(box_relax.molcalc, calc_stress_tensor = true)
        press = get_press(box_relax)
        println("pressure: ", press)
        ## check convergence
        if abs(press) < box_relax.stress_tol
            break
        end
        ## define the direction of volumetric deformation
        new_direc = sign(press*(-1.0))
        ## check if the sign of the stress has changed
        if new_direc != direc
            delta = delta*0.5
            println("delta: ", delta)
        else
            delta = delta*1.01
            println("delta: ", delta)
        end
        direc = new_direc
        ## deform box in that direction
        set_affine_deform_vol!(box_relax.deform, delta*direc, delta*direc, 0.0)
        ## minimize box
        Minimizer = Minimiser(box_relax.molcalc)
        Minimizer.alpha_min = 1.0e0
        Jumol.run_cg!(Minimizer,num_steps=10000, tolerance = box_relax.energy_tol)
    end
end

## biaxial relaxation
function relax_biaxial!(box_relax::Box_relax, delta_x_start::Float64, delta_y_start::Float64; quiet= true)
    if !quiet
        println(" - biaxial relaxation start - ")
    end
    for i in 1:3
        ## relax in x-direction
        relax_direction!(box_relax, 1, delta_x_start, quiet = quiet)
        ## relax in y-direction
        relax_direction!(box_relax, 2, delta_y_start, quiet = quiet)
    end
end

## relaxation in one direction
function relax_direction!(box_relax::Box_relax, num_direc::Int64, delta::Float64; quiet = false)
    ## initial direction of relaxation
    direc = sign(box_relax.molcalc.stress_tensor[num_direc,num_direc]*(-1.0))
    for ii in 1:box_relax.max_num_iter
        if !quiet
            println("relaxation in direction: ", num_direc)
        end
        ## get stress tensor
        Jumol.calc_all_pair_forces!(box_relax.molcalc, calc_stress_tensor = true)
        if !quiet
            println("stress tensor:")
            pretty_table(box_relax.molcalc.stress_tensor,
                        crop = :horizontal, formatters = ft_round(8))
        end
        ## check convergence
        if abs(box_relax.molcalc.stress_tensor[num_direc,num_direc]) < box_relax.stress_tol
            break
        end
        ## define both directions of relaxation
        new_direc = sign(box_relax.molcalc.stress_tensor[num_direc,num_direc]*(-1.0))
        ## check if the sign of the stress has changed
        if new_direc != direc
            delta = delta*0.5
            if !quiet
                println("delta: ", delta)
            end
        else
            delta = delta*1.01
            if !quiet
                println("delta: ", delta)
            end
        end
        direc = new_direc
        ## deform box in that direction
        if num_direc == 1
            set_affine_deform_vol!(box_relax.deform,delta*direc,0.0,0.0)
        end
        if num_direc == 2
            set_affine_deform_vol!(box_relax.deform,0.0,delta*direc,0.0)
        end
        if num_direc == 3
            set_affine_deform_vol!(box_relax.deform,0.0,0.0,delta*direc)
        end
        ## minimize box
        Minimizer = Minimiser(box_relax.molcalc)
        Minimizer.alpha_min = 1.0e0
        Jumol.run_cg!(Minimizer,num_steps=100000,tolerance = box_relax.energy_tol, quiet = quiet)
    end
end

## shear relaxation
function relax_shear!(box_relax::Box_relax, delta::Float64, direc_vec::Vector{Int64})
    ## initial direction of relaxation
    direc = sign(box_relax.molcalc.stress_tensor[direc_vec[1],direc_vec[2]]*(-1.0))
    for ii in 1:box_relax.max_num_iter
        ## get stress tensor
        Jumol.calc_all_pair_forces!(box_relax.molcalc, calc_stress_tensor=true)
        println("relaxation in direction: ", direc_vec)
        println("stress tensor:")
        pretty_table(box_relax.molcalc.stress_tensor,
        	         #noheader = true,
                     crop = :horizontal, formatters = ft_round(8))
        ## check convergence
        if abs(box_relax.molcalc.stress_tensor[direc_vec[1],direc_vec[2]]) < box_relax.stress_tol
            break
        end
        ## define shear direction of relaxation
        new_direc = sign(box_relax.molcalc.stress_tensor[direc_vec[1],direc_vec[2]]*(-1.0))
        ## check if the sign of the stress has changed
        if new_direc != direc
            delta = delta*0.5
            println("delta: ", delta)
        else
            delta = delta*1.01
            println("delta: ", delta)
        end
        direc = new_direc
        ## deform box in that direction
        if direc_vec[1] == 1 && direc_vec[2] == 2
            set_affine_deform_shear!(box_relax.deform,delta*direc, 0., 0.)
        end
        if direc_vec[1] == 1 && direc_vec[2] == 3
            set_affine_deform_shear!(box_relax.deform,0.,delta*direc,0.)
        end
        if direc_vec[1] == 2 && direc_vec[2] == 3
            set_affine_deform_shear!(box_relax.deform,0.,0.,delta*direc)
        end
        ## minimize box
        Minimizer = Minimiser(box_relax.molcalc)
        Minimizer.alpha_min = 1.0e0
        Jumol.run_cg!(Minimizer,num_steps=10000, tolerance = box_relax.energy_tol)
    end
end
