using Distributed
using JLD2
using Plots
using Dates
using BenchmarkTools
using Profile
PBASE = "jumol/Project.toml"
if Base.active_project()[end-length(PBASE)+1:end] != PBASE
    cd(@__DIR__)
    println("Now at ", pwd())
    include("../setup.jl")
    setup_env()
end
using Jumol
Jumol.set_the_current_path(path=@__FILE__)
"""
This demo file demonstrates the aqs procedure which is sped up by a shear in parallel approach
"""

quiet = true
delta_gam = 0.001
max_steps = 100000
eps_mini = 1.0e-6
total_aqs_steps = 400
potential = 33
K = 24
sample_number = 1
temp_struc, molcalc = Jumol.create_struct(sample_number = sample_number)
minier = Jumol.Minimiser(molcalc)

############################
# First the general variables are set and a toy structure with 400 particles is created.
# The shearing will be 0.005 per step, which is the parameter of the transformation Matrix H^{\Delta \gamma}.
# In total the number of AQS steps is fixed to 200 and the minimiser will terminate if the gradient norm is 
# bewlo the eps_mini treshold.
# In the toy structure a binary Lennard Jones potential is used.
############################

@info "NOTE: Depending on your machine this might take a while, since @btime runs the simulation multiple times"
println("We start with the sequentiel procedure")
println("The following timing will show elapsed time and memory usage")
@btime Jumol.aqs_sequentiell!(num_steps = total_aqs_steps, eps=eps_mini, delta_gam = delta_gam, quiet = quiet, minier = deepcopy(minier), potential = potential)
result_dict = Jumol.inner_procedure!(num_steps = total_aqs_steps, eps=eps_mini, delta_gam = delta_gam, quiet = quiet, minier = deepcopy(minier), potential = potential)
println("Final U of sequentiel is ", result_dict["U_hist"][end])
Jumol.plot_U(result_dict["U_hist"], seq=true, only=true)
println()

############################
# Below are extra configurations for the parallel SIP procedure, which is the number of CG steps for the corarse steps
# and the number K of fine steps on the processors as well as the number of processors. 
#IMPORTANT: The following parallel procedure is only faster if the number of threads is larger than 1, 
# for that your Julia session needs to be started with more than one thread, in this case at least 5 are recommenden, 
# paper results were generated with Threads.nthreads()=5.
############################

@info "IMPORTANT: The following parallel procedure is only faster if the number of threads is larger than 1, for that your Julia session needs to be started with more than one thread, in this case at least 5 are recommenden, paper results were generated with Threads.nthreads()=5."
println("The current number of threads is ", Threads.nthreads())
eps_accepter = eps_mini * 3.
P = Threads.nthreads()-1

if P <2
    @error "The number of threads is too low, please start Julia with more threads, e.g. julia -t 5 or set the value in VS Code as julia.NumThreads: auto, this is also added in the .vscode fodler"
else
    minier = Jumol.Minimiser(molcalc, cache_used = true, D_split = K)
    println("Creating shifter object now, and starting parallel procedure") 
    println("The following timing will show elapsed time and memory usage")
    shifter_parallel = Jumol.Aqs_parallel_shifter(deepcopy(minier), potential=potential, #the ibject needs to be availible to all processors using it
                                                            aqs_steps=total_aqs_steps,
                                                            P=P, delta_gam=delta_gam, D_split=K,
                                                            eps_accepter=eps_accepter, 
                                                            max_min_steps=max_steps, eps_mini=eps_mini,
                                                            quiet=quiet)
    println("Starting parallel AQS...")
    @btime U_tracker, final_struct, rts_hist, final_time, overhead_time = Jumol.aqs_sip_procedure_threads(deepcopy(shifter_parallel))
    U_tracker, final_struct, rts_hist, final_time, overhead_time = Jumol.aqs_sip_procedure_threads(shifter_parallel)
    println("Final U of SIP is ", Jumol.final_U(U_tracker, total_aqs_steps))
    Jumol.plot_U(U_tracker, title = "Using K="*string(K), bfig=2600,hfig=1200, only=true)
    println()
end

"""
In the folder samples you can find 1000 uncorelated samples in their respective meta stable state.
"""