## lj force field
# cf. M. Falk, J.S. Langer (1998) PRE 57:7192
################################################
################################################
################################################

mutable struct LJ
    Rc::Float64
    optimize::Bool
    LJ(Rc; optimize = false) = new(Rc, optimize)
    sigma::Vector{Float64}
    epsilon::Vector{Float64}
    V_Rc::Vector{Float64}
    V_prime_Rc::Vector{Float64}
end

## setting all parameters: this must be executed before the potential is used
function set_LJ_params!(lj::LJ; optimize::Bool=false)
    ##
    lj.sigma = MVector{3, Float64}([2*sin(pi/10), 1.0, 2*sin(pi/5)])
    lj.epsilon = MVector{3, Float64}([0.5, 1.0, 0.5])
    #### for shifting and tilding
    ## calculate potential at the cutoff
    lj.V_Rc = zeros(Float64, 3)
    lj.V_Rc[1] = calc_V(lj::LJ, lj.Rc, lj.epsilon[1], lj.sigma[1]);
    lj.V_Rc[2] = calc_V(lj::LJ, lj.Rc, lj.epsilon[2], lj.sigma[2]);
    lj.V_Rc[3] = calc_V(lj::LJ, lj.Rc, lj.epsilon[3], lj.sigma[3]);
    ## calculate derivative of the potential at the cutoff
    lj.V_prime_Rc = zeros(Float64, 3)
    lj.V_prime_Rc[1] = calc_V_prime(lj::LJ, lj.Rc, lj.epsilon[1], lj.sigma[1]);
    lj.V_prime_Rc[2] = calc_V_prime(lj::LJ, lj.Rc, lj.epsilon[2], lj.sigma[2]);
    lj.V_prime_Rc[3] = calc_V_prime(lj::LJ, lj.Rc, lj.epsilon[3], lj.sigma[3]);
    ## optimize speed
    lj.optimize = optimize

end

## calculate the iter-particle potential
function calc_pot(lj::LJ, type1::Int64, type2::Int64, r::Float64)
    if r < lj.Rc
        # current id of atom1 and atom2 -> find parameters
        pos = find_pos(lj,type1,type2)
        sig = lj.sigma[pos]
        eps = lj.epsilon[pos]
        # calculate potential
        return calc_V(lj, r, eps, sig) - lj.V_Rc[pos] - lj.V_prime_Rc[pos]*(r-lj.Rc)
    else
        return 0.0
    end
end

## define pair (find the corresponding positions of the input matrices)
function find_pos(lj::LJ, type_atom1::Int64, type_atom2::Int64)
    # Si-Si
    if type_atom1 == 1 && type_atom2 == 1
        return 1
    # Si-O || O-Si
    elseif (type_atom1 == 1 && type_atom2 == 2) || (type_atom1 == 2 && type_atom2 == 1)
        return 2
    # O-O
    elseif type_atom1 == 2 && type_atom2 == 2
        return 3
    end
    return 1
end

## calculate potential
function calc_V_tester(lj::LJ, dist::Float64, eps::Float64, sig::Float64;u_max::Float64=1.0)
    if dist < sig
        return -u_max/sig * dist + u_max
    else
        return 4 * eps * ((sig/dist)^12.0 - (sig/dist)^6.0)
    end
    #return 4 * eps * ((sig/r)^12.0 ) #
    #return 4 * eps * (-(sig/r)^6.0 ) #
end

function calc_V(lj::LJ, dist::Float64, eps::Float64, sig::Float64)
    return 4 * eps * ((sig/dist)^12.0 - (sig/dist)^6.0)
    #return 4 * eps * ((sig/r)^12.0 ) #
    #return 4 * eps * (-(sig/r)^6.0 ) #
end

## calculate derivative of the potential
function calc_V_prime_tester(lj::LJ, dist::Float64, eps::Float64, sig::Float64;u_max::Float64=1.0)
    if dist<sig
        return - u_max/sig
    else
        return 4 * eps * (6 * sig^6.0/dist^7.0 - 12 *sig^12.0 /dist^13.0)
    end
end

function calc_V_prime(lj::LJ, dist::Float64, eps::Float64, sig::Float64)
    return 4 * eps * (6 * sig^6.0/dist^7.0 - 12 *sig^12.0 /dist^13.0)
end

function calc_V_pp(dist::Float64, eps::Float64, sig::Float64)
    return 4 * eps * (- 42 * sig^6.0/dist^8.0 + 156 *sig^12.0 /dist^14.0)
end