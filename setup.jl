using Pkg

function setup_env()
    println("Setting up environment")
    if pwd()[end-5:end] != Base.active_project()[1:end-13]
        println("Need to change to project root in order to be on correct level to Project.toml")
        cd(Base.active_project()[1:end-13])
        println("Now at ", pwd())
    end
    Pkg.instantiate();
    Pkg.precompile();
    Pkg.activate(".");
    Pkg.resolve();
    Pkg.instantiate(); 
    #Pkg.test();
end