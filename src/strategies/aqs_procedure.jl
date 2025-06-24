"""
This gives an aqs procedure object
currrently only supports simple shearing in one direction and for two dimensions
"""

#TODO make availible for 3D
"""
    Aqs_shifter
    This object defines a possible shift of given strucure, it only needs a calc object. All options on how to shift and do the aqs procdure is defined in the object.
"""
mutable struct Aqs_shifter
    calc::Calc
    deformer::Aff_deform
    minimiser::Minimiser
    aqs_steps::Int64
    min_method::String
    alpha_min::Float64
    acc_factor::Float64
    max_min_steps::Int64
    eps::Float64
    delta_gam::Float64
    dist_break::Float64
    linear_guessing::Bool
    q_num::Int64
    reg::Float64
    quiet::Bool
    visualize::Bool
    siO_glass::Bool

    function Aqs_shifter(calc;
        aqs_steps = 30,
        min_method= "cg",
        alpha_min = 1.0e0,
        acc_factor = 1.001,
        max_min_steps=100000,
        eps = 1.0e-5,
        delta_gam = 0.005,
        dist_break = 0.5,
        linear_guessing= true,
        q_num=4,
        reg=0.0,
        quiet = true,
        visualize = false,
        siO_glass = false)
        minimiser = Jumol.Minimiser(calc)
        deformer = Jumol.Aff_deform(calc.structure)
        minimiser.alpha_min = alpha_min
        minimiser.acc_factor = acc_factor
        if siO_glass
            delta_gam = delta_gam*100.
            dist_break = dist_break * 2.
        end
        return new(calc, deformer, minimiser, 
        aqs_steps, min_method, alpha_min, acc_factor, max_min_steps, eps, delta_gam, dist_break, linear_guessing, q_num, reg, quiet, visualize, siO_glass)
    end
end

"""
    U_obj
    This object gives the opportunity to track all the solutions calculated and if they are accepted or rejected. There are also informational values helping to
    identify and analyse results.
"""
mutable struct U_obj
    U::Float64
    two_norm::Float64
    relative_index::Int64
    coarse_index::Int64
    abs_index::Int64
    given_link::Tuple{Int64, Int64}
    steps::Int64
    stress::Float64
    accepted::Bool
    U_obj(U, two_norm, relative_index, coarse_index, abs_index, given_link, steps, stress; accepted = false) = new(U, two_norm, relative_index, coarse_index,abs_index, given_link, steps, stress, accepted)
end

"""
    get_positions_vector(structure::Structure)
    Returns a vector of length 3*structure.noa, note that all information about shear of the box, stress and forces is lost if just the positions are saved
"""
function get_positions_vector(structure::Structure, posis::Vector{Float64})
    for i in eachindex(structure.atom_list)
        posis[(i-1)*3+1:i*3] .= structure.atom_list[i].pos
    end
end

"""
    push_matrix!(mat::Matrix{Float64}, vec::Matrix{Float64})
    Pushes vec to the front of given matrix, length of vector must be number of rows of matrix 
"""
function push_matrix!(mat::Matrix{Float64}, vec::Vector{Float64})
    cs = size(mat,2)
    @views mat[:, 1:end-1] .= mat[:,2:end]
    @views mat[:,end] .= vec
end

"""
    trafo_and_normalize(A::Matrix{Float64}, B::Matrix{Float64})
    Centers and normalizes the matrix in order to make the alpha vector which is calculated to make a linear initial guess more interpretable.
"""
function trafo_and_normalize(A::Matrix{Float64}, B::Matrix{Float64})
    temp_A = (A .- mean(A))
    temp_B = (B .- mean(A)) # trafo by only the mean of A because need same trafo 
    return temp_A ./ maximum(abs.(temp_A)), temp_B ./ maximum(abs.(temp_A)), mean(A), maximum(abs.(temp_A))
end

"""
    trafo_back(V::Vector{Float64}, shift::Float64, scale::Float64)
    Reverses the trafo_and_normalize to shift back to the original scale
"""
function trafo_back!(V::Vector{Float64}, shift::Float64, scale::Float64)
    V.= (V.*scale) .+ shift 
end

"""
    is_plastic_event(rts::Matrix{Float64}; tolerance = 10.0)
    Checks in the history if there was a plastic event, because then linear guessing does not make any sense and we need to start again building the history and not do any linear guessing for the next q steps.
        Check is done only by comparing subsequent norms to a given tolerance.
"""
function is_plastic_event(rts::AbstractArray{Float64}; tolerance = 10.0)
    for k in 1:size(rts,2)-2
        @views fac = norm(rts[:,k]-rts[:,k+1])/norm(rts[:,k+1]-rts[:,k+2])
        if fac > tolerance
            return true
        end
    end
    return false
end

"""
    construct_init_guess!(r_before::Matrix{Float64}, r_after::Matrix{Float64}; reg = 0., exp_reg = nothing, two_dim = false)
    Given the hisroty of rt and rts we regress linear factors and use the first vector of rt (rest of matrix is not needed but easier to keep track of rt). Method only works if there is no incontinouity before.
    This is checked by is_plastic_event. A new inital guess for the rt is returend which cosiderably lowers the number of minisation steps and hence speeds up the procedure.
"""
function construct_init_guess!(r_before::AbstractArray{Float64}, r_after::AbstractArray{Float64}, r_predict::Vector{Float64}; ent_p_relative::Int64=0)
    # as opposed to the notation in the paper r_before is the configs after affine transformation but before minimization
    # meaning r_after is the config history after minimization
    if is_plastic_event(r_before)
        @views r_predict .= r_before[:,end]
    else
        @views begin
            r_temp = (r_before[:, 1+ent_p_relative:end-1] .- r_after[:,2+ent_p_relative])
            alpha = (transpose(r_temp) * r_temp)\(transpose(r_temp) * (r_before[:,end] - r_after[:,2+ent_p_relative]))
            mul!(r_predict, r_after[:,2+ent_p_relative:end], alpha)
            r_predict .+= (1-sum(alpha)) .* r_after[:,2+ent_p_relative]
            if norm(abs.(r_predict-r_before[:,end]))>0.1
                @views r_predict .= r_before[:,end]
            end
        end
    end
end

function basis_poly(delta_index::Int64, deltas::Vector{Float64})
    L_i = delta::Float64 -> prod((delta-delta_j)/(deltas[delta_index]-delta_j) for delta_j in deltas if delta_j !=deltas[delta_index])
    return L_i
end

function lagrange_extrapolation!(extrapol::Vector{Float64}, r_before::AbstractArray{Float64}, r_after::AbstractArray{Float64}, total_gam::Float64, delta_gam::Float64)
    if is_plastic_event(r_before)
        @views extrapol .= r_before[:,end]
    else
        N3, q_num = size(r_after)
        deltas = reverse(Float64[total_gam - delta_gam * i for i in 0:q_num-1])
        extrapol .= zeros(Float64, N3)
        extrapolation_gam = total_gam + delta_gam
        for q in 1:q_num 
            extrapol .+= r_after[:,q] * basis_poly(q ,deltas)(extrapolation_gam)
        end
        if norm(abs.(extrapol-r_before[:,end]))>0.1
            @views extrapol .= r_before[:,end]
        end
    end
end

function hesse_extrapolation!(init_guess::Vector{Float64}, r_before::AbstractArray{Float64}, basis::Matrix{Float64}, alpha::Vector{Float64})
    if is_plastic_event(r_before)
        @views init_guess .= r_before[:,end]
        return false
    else
        @views init_guess .= r_before[:,end] + basis * alpha
        if norm(abs.(init_guess-r_before[:,end]))>0.1
            @views init_guess .= r_before[:,end]
            return false
        end
        return true
    end
end

function del_zeros(veci :: Vector{Float64})
    temp_no = zeros(Float64, Int64(2/3 * length(veci)))
    noa = Int64(length(veci)/3)
    for i in 1:noa
        temp_no[2*i-1:2*i] = veci[3*i-2:3*i-1]
    end
    return temp_no
end

function del_zeros(veci :: Matrix{Float64})
    temp_no = zeros(Float64, Int64(2/3 * size(veci,1)), size(veci,2))
    noa = Int64(size(veci,1)/3)
    for i in 1:noa
        temp_no[2*i-1:2*i, :] = veci[3*i-2:3*i-1, :]
    end
    return temp_no
end

function add_zeros(veci :: Vector{Float64})
    temp_no = zeros(Float64, Int64(3/2* length(veci)))
    noa = Int64(length(veci)/2)
    for i in 1:noa
        temp_no[3*i-2:3*i-1] .= veci[2*i-1:2*i]
        temp_no[3*i] = 0.0
    end
    return temp_no
end


"""
    set_positions!(structure, positions::Vector{Float64})
    Sets position given a structure and positions, does change the structure but does NOT change the geometry accordingly. ONLY use for plotting purposes.
"""
function set_positions!(structure::Jumol.Structure, positions::Vector{Float64})
    """
    setting postions in order to plot and visualize better
    NOTE: This is only for visualization and changes the struct, the potential energy and forces are not updated correctly here, also the box might be missaligned if the inital strucutre was sheared further
    """;
    for n in 1:structure.noa
        structure.atom_list[n].pos .= positions[3*(n-1)+1:3*n]
    end
end

"""
    plot_net!(structure; posis = nothing, title_text = nothing, shear_step = nothing, size1 = 12.0, size2 = 22.5)
    Plots a netwrok given size, can change positions in order to show AQS progress, will change strucutre, need deeepcopy if want to use later.
"""
function plot_net!(structure::Jumol.Structure; posis::Union{Nothing, Vector{Float64}} = nothing, title_text = nothing, shear_step = nothing, size1::Float64=12.0, size2::Float64=22.5, fig_size::Number=550, paper::Bool=false)
    """
    visualizes a net given, can give position vector which changes the strucuture, number of steps can be calculated from the box shear of the structure
    """;
    if !isnothing(posis)
        set_positions!(structure, posis)
    end

    hfig = fig_size
    bfig = fig_size
    Plotter = Jumol.Vis2d(structure)
    Plotter.bfig = bfig
    Plotter.hfig = hfig
    fig = Jumol.plot_box(Plotter)
    Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                        [-structure.box.lx*0.075, structure.box.lx*1.40, structure.box.lx*1.40,-structure.box.lx*0.075],
                        [-structure.box.ly*0.075, -structure.box.ly*0.075, structure.box.ly*1.075, structure.box.ly*1.075] )
    if !isnothing(title_text)
        title!(title_text)
    elseif !isnothing(shear_step) 
        step = Int64(round(structure.box.lxy / shear_step))
        title!("AQS step "*string(step))
    end
    
    if paper
        quiver!([0], [0], quiver=([structure.box.lx],[0.0]), color="black",lw=3)
        quiver!([0], [0], quiver=([0.0], [structure.box.ly]), color="black",lw=3)

        # Add annotation text near the arrow
        annotate!(-3, structure.box.ly-3, text(L"$\mathbf{\hat{h}_1}$", :black, 20))
        annotate!(structure.box.lx-3, -3, text(L"$\mathbf{\hat{h}_2}$", :black, 20))
    end

    display(fig)
    return fig
end

"""
    push_dict_to_U!(U_tracker, info_dict, absolut_tracker, fine_steps, boarder)
    Push all information from Dicts to U_objects with corresponding statuses to keep track of all solutions.
"""
function push_dict_to_U!(U_tracker::Vector{U_obj}, info_dict::Dict{Int64,Dict{String, Union{Array{Union{Float64, Int64}}, Nothing}}}, absolut_tracker::Int64, fine_steps::Int64, boarder::Int64)
    """
    gets data of all parallel aqs steps from the info dict and creates a u_obj for every aqs step and marks them as accepted or denied
    """;
    for coarse_k in keys(info_dict)
        if coarse_k<=boarder
            status = true
        else
            status = false
        end

        giver = (coarse_k-1, absolut_tracker) #-1 in order to make it comform with plots
        for k in 1:length(info_dict[coarse_k]["U"])
            if coarse_k > 1
                temp_absolut_tracker = absolut_tracker - 1
            else
                temp_absolut_tracker = absolut_tracker
            end
            #println("This is ",k, " ", coarse_k, " ", giver, " ", status)
            @views begin
                push!(U_tracker, U_obj(info_dict[coarse_k]["U"][k], 
                                norm(info_dict[coarse_k]["rts"][:,k]),
                                k, 
                                coarse_k,
                                temp_absolut_tracker + fine_steps*(coarse_k-1)+ k,
                                giver,
                                info_dict[coarse_k]["steps"][k],
                                info_dict[coarse_k]["stress"][k],
                                accepted = status))
            end
        end
    end
end

"""
    plot_U(U_tracker; bfig = 2200, hfig = 1200, scale = 50., seq = false)
    Plots the solutions which are rejected and accepted in correspondign color, can also plot the sequentiell solutions.
"""
function plot_U(U_tracker::Union{Vector{U_obj},Vector{Float64}}; bfig::Int64= 2200, hfig::Int64 = 1200, scale::Float64 = 50., seq::Bool=false, title=nothing, stacked::Bool=false, delta_gam::Float64 = 0.005, ly::Float64 = 50.0, only::Bool=false)
    """
    plots all u_obj to a graph relativ to the aqs step and indicates first steps after coarse, and also if accepted or denied
    """;
    fsize=55
    marker_orange = 25
    fsize_marker = 22
    marker_green = 8
    marker_red = 8

    if seq 
        points = []
        for k in eachindex(U_tracker)
            point = (k, U_tracker[k])
            push!(points, point)
        end
        fsize=55
        yformatter = x -> @sprintf("%.3e", x)
        x_vals = [points[k][1]*delta_gam/ly for k in eachindex(points)]
        y_vals = [points[k][2] for k in eachindex(points)]
        y_n = minimum(y_vals)
        y_x = maximum(y_vals)
        yticks_custom = range(y_n + (y_x-y_n)*0.2, stop= y_x - (y_x-y_n)*0.1, length=3)
        if only
            bmargin = 25Plots.mm
        else
            bmargin = -15Plots.mm
        end
        p=plot(x_vals, y_vals,legend=false, xticks=false, xgrid=true, xlabel=!stacked, color = "black", lw=9, size=(bfig, hfig), label = "Accepted steps",margin=margin=15Plots.mm, xtickfontsize = fsize, ytickfontsize =fsize-20, legendfontsize =fsize, titielfontsize= fsize, xlabelfontsize=fsize, ylabelfontsize=fsize, titlefontsize=fsize, framestyle=:box, top_margin = 15Plots.mm, bottom_margin=bmargin, yformatter = yformatter, yticks = yticks_custom)
        ylabel!(L"$\mathbf{\mathcal{U}(r,\gamma)}$")
        if !stacked
            xlabel!("AQS Steps")
        end
        if !isnothing(title)
            title!(title)
        end
        display(p)
        return p
    else
        medi = median([val.U for val in U_tracker])
        after_coarse = []
        accepted = []
        rejected = []
        for val in U_tracker
            if medi - scale < val.U < medi + scale
                point = (val.abs_index, val.U, val.coarse_index -1 )
                if val.relative_index == 1 
                    push!(after_coarse, point)
                elseif val.accepted 
                    push!(accepted, point)
                else
                    push!(rejected, point)
                end
            end
        end
        yformatter = x -> @sprintf("%.3e", x)
        x_vals = [accepted[k][1]*delta_gam/ly for k in eachindex(accepted)]
        x_rej = [rejected[k][1]*delta_gam/ly for k in eachindex(rejected)]
        x_ac = [after_coarse[k][1]*delta_gam/ly for k in eachindex(after_coarse)]
        y_vals = [accepted[k][2] for k in eachindex(accepted)]
        y_n = minimum(y_vals)
        y_x = maximum(y_vals)
        yticks_custom = range(y_n + (y_x-y_n)*0.2, stop= y_x - (y_x-y_n)*0.1, length=3)
        if only
            bmargin = 25Plots.mm
        else
            bmargin = -15Plots.mm
        end
        p=scatter(x_vals, y_vals, xticks=!stacked, xlabel=!stacked, color = "green", ms = marker_green, size=(bfig, hfig), label = "Accepted steps",margin=margin=15Plots.mm, xtickfontsize = fsize, ytickfontsize =fsize-20, legendfontsize =fsize-25, titielfontsize= fsize, xlabelfontsize=fsize, ylabelfontsize=fsize, titlefontsize=fsize, framestyle=:box, top_margin = 15Plots.mm, bottom_margin=bmargin, yformatter = yformatter, yticks = yticks_custom)
        scatter!(x_rej, [rejected[k][2] for k in eachindex(rejected)], color = "red", ms = marker_red, size=(bfig, hfig), label = "Rejected steps")
        scatter!(x_ac, [after_coarse[k][2] for k in eachindex(after_coarse)], color = "orange", ms = marker_orange, size=(bfig, hfig), label = "First step after coarse guess",
        series_annotations = text.([after_coarse[k][3] for k in eachindex(after_coarse)], fsize_marker))
        ylabel!(L"${\mathcal{U}(\mathbf{r},\gamma)}$")
        if !stacked
            xlabel!("AQS Steps")
        end
        
        if !isnothing(title)
            title!(title)
        end
        display(p)
        return p
    end
end

"""
    plot_stress(U_tracker; bfig = 2200, hfig = 1200, scale = 50.)
    Plots the stress over the course of the AQS procedure. Similar to plot_U.
"""
function plot_stress(U_tracker::Union{Vector{U_obj},Vector{Float64}}; bfig::Int64= 2200, hfig::Int64 = 1200, scale::Float64 = 50., legend::Bool=true, seq::Bool=false, delta_gam::Float64 = 0.005, ly::Float64 = 50.0)
    """
    same as above but plots the average stress
    """;
    fsize=55
    marker_orange = 25
    fsize_marker = 22
    marker_green = 8
    marker_red = 8
    if seq 
        points = []
        for k in eachindex(U_tracker)
            point = (k, U_tracker[k])
            push!(points, point)
        end
        fsize=55
        yformatter = x -> @sprintf("%.3e", x)
        x_vals = [points[k][1]*delta_gam/ly for k in eachindex(points)]
        y_vals = [points[k][2] for k in eachindex(points)]
        y_n = minimum(y_vals)
        y_x = maximum(y_vals)
        yticks_custom = range(y_n + (y_x-y_n)*0.1, stop=y_x - (y_x-y_n)*0.2, length=3)
        p=plot(x_vals, y_vals, xticks=true, legend=false, color = "black", lw=9, size=(bfig, hfig), margin=15Plots.mm, xtickfontsize = fsize-20, ytickfontsize =fsize-20, legendfontsize =fsize, titielfontsize= fsize, xlabelfontsize=fsize, ylabelfontsize=fsize, titlefontsize=fsize, framestyle=:box, top_margin = -12Plots.mm, yformatter = yformatter, yticks = yticks_custom)
        ylabel!(L"$\mathbf{\tau}$")
        xlabel!(L"$\mathbf{\gamma}$")
        
        display(p)
        return p
    else
        medi = median([val.stress for val in U_tracker])
        after_coarse = []
        accepted = []
        rejected = []
        for val in U_tracker
            if medi - scale < val.stress < medi + scale
                point = (val.abs_index, val.stress, val.coarse_index -1 )
                if val.relative_index == 1 
                    push!(after_coarse, point)
                elseif val.accepted 
                    push!(accepted, point)
                else
                    push!(rejected, point)
                end
            end
        end
        yformatter = x -> @sprintf("%.3e", x)
        x_vals = [accepted[k][1]*delta_gam/ly for k in eachindex(accepted)]
        x_rej = [rejected[k][1]*delta_gam/ly for k in eachindex(rejected)]
        x_ac = [after_coarse[k][1]*delta_gam/ly for k in eachindex(after_coarse)]
        y_vals = [accepted[k][2] for k in eachindex(accepted)]
        y_n = minimum(y_vals)
        y_x = maximum(y_vals)
        yticks_custom = range(y_n + (y_x-y_n)*0.1, stop=y_x - (y_x-y_n)*0.2, length=3)
        p=scatter(x_vals, y_vals, color = "green", ms = marker_green, size=(bfig, hfig), label = "Accepted steps",margin=margin=15Plots.mm, xtickfontsize = fsize-20, ytickfontsize =fsize-20, legendfontsize =fsize, titielfontsize= fsize, xlabelfontsize=fsize, ylabelfontsize=fsize, titlefontsize=fsize, framestyle=:box, legend=legend, yformatter = yformatter, yticks = yticks_custom,  top_margin = -9Plots.mm,)
        scatter!(x_rej, [rejected[k][2] for k in eachindex(rejected)], color = "red", ms = marker_red, size=(bfig, hfig), label = "Rejected steps")
        scatter!(x_ac, [after_coarse[k][2] for k in eachindex(after_coarse)], color = "orange", ms = marker_orange, size=(bfig, hfig), label = "First step after coarse guess",
        series_annotations = text.([after_coarse[k][3] for k in eachindex(after_coarse)], fsize_marker))
        ylabel!(L"$\tau$")
        xlabel!(L"\gamma")
        display(p)
        return p
    end
end

function plot_U_and_stress(U_tracker::Union{Vector{U_obj},Vector{Vector{Float64}}}; bfig::Int64= 2200, hfig::Int64 = 1200, scale::Float64 = 50., seq::Bool=false, title=nothing)
    """
    plots all u_obj to a graph relativ to the aqs step and indicates first steps after coarse, and also if accepted or denied
    """;
    if seq
        p_U = plot_U(U_tracker[1], bfig = bfig, hfig =hfig, scale = scale, stacked=true, title=title, seq = seq)
        p_stress = plot_stress(U_tracker[2], bfig = bfig, hfig =  hfig, scale = scale, legend = false, seq = seq)
        p_all = plot(p_U, p_stress, layout = @layout([a; b]), link = :x, margin=25Plots.mm)
    else
        p_U = plot_U(U_tracker, bfig = bfig, hfig =hfig, scale = scale, stacked=true, title=title, seq = seq)
        p_stress = plot_stress(U_tracker, bfig = bfig, hfig =  hfig, scale = scale, legend = false, seq = seq)
        p_all = plot(p_U, p_stress, layout = @layout([a; b]), link = :x)
    end

    display(p_all)
    return p_all
end

"""
    final_U(U_tracker, total_steps)
    Gets final U_obj given total steps. Maybe there is to many solutions hence we need to cut the ones which are to long.
"""
function final_U(U_tracker::Vector{U_obj}, total_steps::Int64)
    """
    returns the final U value which was accepted and is exactly the number of steps
    """;
    index = [k for k in U_tracker if k.abs_index == total_steps&&k.accepted]
    if length(index) == 0
        println("No value given index")
        return 0.
    else
        return index[1].U
    end
end

"""
    zoom_plot!(x_range, y_range)
    Function to zoom into plots in order to see details better.
"""
function zoom_plot!(x_range::Tuple{Int64,Int64}, y_range::Tuple{Int64,Int64})
    """
    zooms into plot in given range
    """;
    xlims!(x_range)
    ylims!(y_range)
end

"""
    plot_mini(min_hist;level = 0, bfig = 1100, hfig = 600)
    Plots the minimisation in a singular AQS step, can plot more than one minimisation into one plot.
"""
function plot_mini(min_hist::Vector ;level::Int64 = 0, bfig::Int64 = 1100, hfig::Int64 = 600)
    """
    plots minimisation history given a certain step, mostly used for coarse adaptive search to interepret results
    """;
    p = plot(1:length(min_hist[1]), min_hist[1], label = "1. aqs step", size=(bfig, hfig), margin=margin=10Plots.mm)
    for k in 2:length(min_hist)
        plot!(1:length(min_hist[k]), min_hist[k], label = string(k)*". aqs step")
    end
    #vline!([coarse_max_cg], color = "red", label = "Max coarse cg steps")
    title!("Minimisation of coarse steps on level "*string(level))
    xlabel!("CG steps")
    ylabel!(L"Potential Energy $\mathbf{\mathcal{U}}$")
    display(p)
end


"""
    plot_solution(solution_dict::Dict)
    This function plots all the coarse steps.
"""
function plot_solution_temp(solution_dict::Dict; title_text = nothing)
    temp_U_cache = []
    for ind in eachindex(solution_dict["U_hist"])
        temp_obj = U_obj(solution_dict["U_hist"][ind],
            0.,
            ind,
            0,
            ind,
            (0,0),
            0,
            0; accepted = true)
        push!(temp_U_cache, temp_obj)
    end
    plot_U(temp_U_cache, title = title_text)
end

function plot_U_gif(U_tracker::Vector{Jumol.U_obj}; bfig = 2200, hfig = 1200, title=nothing, fsize=25)
    P_num = maximum([k.coarse_index for k in U_tracker])
    D_split = maximum([k.relative_index for k in U_tracker])
    num_steps = maximum([k.abs_index for k in U_tracker])
    u_min = minimum([k.U for k in U_tracker]) - 2.0
    u_max = maximum([k.U for k in U_tracker]) + 2.0
    starters= sort(unique([k.given_link[2] for k in U_tracker]))
    relative_counter = -1
    coarse_counter = -1
    still_plotting = true
    coarse_this_done = false
    coarse_relative = 1
    fine_relative = 2
    anim = @animate while still_plotting
        if coarse_counter == -1
            p = plot([],[], xlims = (0,num_steps), ylims=(u_min, u_max), primary=false, size=(bfig, hfig),margin=margin=15Plots.mm, xtickfontsize = fsize, ytickfontsize =fsize, legendfontsize =fsize, titielfontsize= fsize, xlabelfontsize=fsize, ylabelfontsize=fsize, titlefontsize=fsize)
            ylabel!(L"Potential Energy $\mathbf{\mathcal{U}}$")
            xlabel!("AQS Steps")
            if !isnothing(title)
                title!(title)
            end
            coarse_counter = 1
        else 
            all_points = [k for k in U_tracker if starters[coarse_counter] == k.given_link[2]]
            if !coarse_this_done 
                x_points_coarse = [k.abs_index for k in all_points if k.relative_index==1&&k.coarse_index==coarse_relative]
                y_points_coarse = [k.U for k in all_points if k.relative_index==1&&k.coarse_index==coarse_relative]
                scatter!(x_points_coarse, y_points_coarse, ms = 15, color = "orange", label="Coarse starting at "*string(starters[coarse_counter]),primary=false)
                annotate!(num_steps, u_max, text("█"^11, :white, :right, 60))
                annotate!(num_steps, u_max, text("█"^11, :white, :left, 60))
                annotate!(num_steps, u_max, text("Sequentiell", :orange, :right, 40))
                if coarse_relative == P_num
                    coarse_this_done = true
                    coarse_relative = 1
                else
                    coarse_relative += 1
                end
            else
                x_inner = [k.abs_index for k in all_points if k.relative_index==fine_relative]
                y_inner = [k.U for k in all_points if k.relative_index==fine_relative]
                scatter!(x_inner, y_inner, ms = 10, color = "black", primary=false)
                annotate!(num_steps, u_max, text("█"^11, :white, :right, 60))
                annotate!(num_steps, u_max, text("█"^11, :white, :left, 60))
                annotate!(num_steps, u_max, text("Parallel", :green, :right, 40))
                if fine_relative == D_split
                    all_points = [k for k in U_tracker if starters[coarse_counter] == k.given_link[2]]
                    x_acc = [k.abs_index for k in all_points if k.accepted==true]
                    y_acc = [k.U for k in all_points if k.accepted==true]
                    x_rej = [k.abs_index for k in all_points if k.accepted==false]
                    y_rej = [k.U for k in all_points if k.accepted==false]
                    scatter!(x_acc, y_acc, ms = 10, color = "green", primary=false)
                    scatter!(x_rej, y_rej, ms = 10, color = "red", primary=false)
                    annotate!(num_steps, u_max, text("■"^11, :white, :right, 40))
                    annotate!(num_steps, u_max, text("■"^11, :white, :left, 40))
                    annotate!(num_steps, u_max, text("Parallel", :green, :right, 40))
                    fine_relative = 2
                    coarse_this_done = false
                    if coarse_counter >= length(starters)
                        still_plotting = false
                    else
                        coarse_counter += 1
                    end
                else
                    fine_relative += 1
                end
            end

        end
    end

    gif(anim, "analysis/procedure_example.gif", fps = 2)
end