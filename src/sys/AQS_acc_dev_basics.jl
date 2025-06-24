
flush(stdout)

function build_basis(pairs::Vector{Tuple{Float64, Vector{Float64}}}, k::Int64; to_3d::Bool=true)
    if k<=2
        return zeros(0,0)
    else
        if to_3d
            size_3d = Int64(3*size(pairs[1][2],1)/2)
            base_ks = zeros(size_3d, k-2)
            for i in 1:k-2
                for k in 1:Int64(size_3d/3)
                    base_ks[3*k-2:3*k-1, i] .= pairs[i+2][2][2*k-1:2*k]
                    base_ks[3*k, i] = 0.0
                end
            end
        else
            base_ks = zeros(size(pairs[1][2],1), k-2)
            for i in 1:k-2
                base_ks[:,i] .= pairs[i+2][2]
            end
        end

    return base_ks
    end
end 

function build_lmbd_vec(pairs::Vector{Tuple{Float64, Vector{Float64}}}, k::Int64)
    lmbd_vec = Float64[ep[1] for ep in pairs[3:k]]
    return lmbd_vec
end

function fit_alphas(ev::Vector{Float64}, base::Matrix{Float64})

    A = transpose(base) * base
    b = transpose(base) * ev
    alpha = A\b
    build_ev = base * alpha

    return alpha, build_ev
end

function basis_poly(delta_index::Int64, deltas::Vector{Float64})
    L_i = delta::Float64 -> prod((delta-delta_j)/(deltas[delta_index]-delta_j) for delta_j in deltas if delta_j !=deltas[delta_index])
    return L_i
end


function create_struct(;sample_number::Int64=1)
    """
    create a toy stuct to test methods, a LJ glass and also SIO glass for more complex testing
    """;
    molstruc = Structure()
    molstruc.rc = 2.4
    molstruc.rskin = 0.4
    molstruc.pbx = 1
    molstruc.pby = 1
    molstruc.linked_cells_bool = true
    
    ## load LJ glass
    file_name = "../samples/cooled_samples/sample_cooled_20_20_no_"*string(sample_number)*".lammpstrj"
    Jumol.read_lammpstrj!(molstruc,file_name)
    Jumol.initialize_structure!(molstruc)
    molcalc = Jumol.Calc(molstruc)
    Jumol.initialize_potential!(molcalc,33)


    return molstruc, molcalc
end

