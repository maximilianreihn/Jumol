## Yukawa force field (Roy, 2018)
######################################
######################################
######################################

mutable struct Yukawa_2D
    Rc::Float64
    kappa::Float64
    Yukawa_2D(Rc; kappa = 1.0) = new(Rc, kappa)
    params::Vector{Vector{Float64}}
    V_Rc::Vector{Float64}
    V_prime_Rc::Vector{Float64}
end

## setting all parameters: this must be executed before the potential is used
function set_Yukawa_2D_params!(yukawa::Yukawa_2D)
    #### parameter setting according to Roy, 2018
    ## row 1: sigma
    ## row 2: q
    ## col 1: Si-Si , col 2: Si-O, col 3: O-O
    #TODO map/dict
    yukawa.params = Vector{Vector{Float64}}([Vector{Float64}([2.250,1.075,0.900]), Vector{Float64}([1.500,-1.000,0.670])])
    #### long range parameter kappa
    ## Carre et al. (2007) The journal of chemical physics 127:114512.
    yukawa.kappa = 1.0/5.649
    #### for shifting and tilding
    ## calculate potential at the cutoff
    yukawa.V_Rc = zeros(Float64, 3)
    yukawa.V_Rc[1] = calc_V(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1][1],yukawa.params[2][1]); #sigma(Si-Si), q(Si-Si)
    yukawa.V_Rc[2] = calc_V(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1][2],yukawa.params[2][2]); #Sigma(Si-O) , q(Si-O)
    yukawa.V_Rc[3] = calc_V(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1][3],yukawa.params[2][3]); #Sigma(O-O)  , q(O-O)
    ## calculate derivative of the potential at the cutoff
    yukawa.V_prime_Rc = zeros(Float64, 3)
    yukawa.V_prime_Rc[1] = calc_V_prime(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1][1],yukawa.params[2][1]); #sigma(Si-Si), q(Si-Si)
    yukawa.V_prime_Rc[2] = calc_V_prime(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1][2],yukawa.params[2][2]); #Sigma(Si-O) , q(Si-O)
    yukawa.V_prime_Rc[3] = calc_V_prime(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1][3],yukawa.params[2][3]); #Sigma(O-O)  , q(O-O)
end

## calculate the inter-particle potential
function calc_pot(yukawa::Yukawa_2D, type1::Int64, type2::Int64, dist::Float64)
    pot = 0.0
    if dist < yukawa.Rc
        # current id of atom1 and atom2 -> find parameters
        pos = find_pos(yukawa,type1,type2)
        # calculate potential
        p1 = yukawa.params[1][pos]
        p2 = yukawa.params[2][pos]
        vrc = yukawa.V_Rc[pos]
        vrcp = yukawa.V_prime_Rc[pos]
        @views pot = calc_V(yukawa, dist, p1,p2) - vrc - vrcp*(dist-yukawa.Rc)
    end
  return pot
end

## define pair (find the corresponding positions of the input matrices)
function find_pos(yukawa::Yukawa_2D, type_atom1::Int64, type_atom2::Int64)
    col = 1
    # Si-Si
    if type_atom1 == 1 && type_atom2 == 1
        col = 1
    end
    # Si-O
    if type_atom1 == 1 && type_atom2 == 2
        col = 2
    end
    # O-Si
    if type_atom1 == 2 && type_atom2 == 1
        col = 2
    end
    # O-O
    if type_atom1 == 2 && type_atom2 == 2
        col = 3
    end
    return col
end


function calc_V_pp(yukawa::Yukawa_2D, r::Float64, sig::Float64, q::Float64)
    return 1.0
end

## calculate potential (optimized version)
function calc_V(yukawa::Yukawa_2D, r::Float64, sig::Float64, q::Float64)
    one_over_r = 1.0/r
    sig_over_r = one_over_r*sig
    sig_over_r2 = sig_over_r*sig_over_r
    sig_over_r4 = sig_over_r2*sig_over_r2
    sig_over_r8 = sig_over_r4*sig_over_r4
    q_over_r = one_over_r*q
    return sig_over_r8*sig_over_r4 + q_over_r*exp(-yukawa.kappa*r)
end

## calculate derivative of the potential (optimized version)
function calc_V_prime(yukawa::Yukawa_2D, r::Float64, sig::Float64, q::Float64)
    one_over_r = 1/r
    q_over_r = one_over_r*q
    one_over_r2 = one_over_r*one_over_r
    sig_over_r2 = one_over_r2*sig*sig
    sig_over_r4 = sig_over_r2*sig_over_r2
    sig_over_r8 = sig_over_r4*sig_over_r4
    q_over_r2 = one_over_r2*q
    return (-12.0)*sig_over_r8*sig_over_r4*one_over_r - yukawa.kappa*q_over_r * exp(-yukawa.kappa*r) - q_over_r2 * exp(-yukawa.kappa*r)
end

## Yukawa force field (Roy, 2018)
######################################
######################################
######################################
