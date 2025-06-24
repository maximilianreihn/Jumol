## system functions
######################################
######################################
######################################
#
# (c) Franz Bamer, Mai-2022
######################################


## set the current directory for calculations
function set_the_current_path(;path::String=Base.source_path())
	path = collect(path)
	while path[end] != '/'
		char_del = path[end]
		num_delete = length(path)
		path = deleteat!(path,num_delete)
	end
	path_to_current_file = String(path)
	cd(path_to_current_file)
	#println("switching to")
	run(`pwd`)
end
