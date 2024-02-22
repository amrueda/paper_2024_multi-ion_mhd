#Script to run convergence tests

# Convergence test with EC fluxes
open(joinpath(@__DIR__,"log_ec_n2.txt"), "w") do io
     redirect_stdout(io) do
         convergence_test(joinpath(@__DIR__,"elixir_convergence_ec.jl"), 4, polydeg=2, initial_refinement_level = 4)
     end
end
open(joinpath(@__DIR__,"log_ec_n3.txt"), "w") do io
     redirect_stdout(io) do
         convergence_test(joinpath(@__DIR__,"elixir_convergence_ec.jl"), 4, polydeg=3, initial_refinement_level = 3)
     end
end
open(joinpath(@__DIR__,"log_ec_n4.txt"), "w") do io
     redirect_stdout(io) do
         convergence_test(joinpath(@__DIR__,"elixir_convergence_ec.jl"), 4, polydeg=4, initial_refinement_level = 3)
     end
end
open(joinpath(@__DIR__,"log_ec_n5.txt"), "w") do io
     redirect_stdout(io) do
         convergence_test(joinpath(@__DIR__,"elixir_convergence_ec.jl"), 4, polydeg=5, initial_refinement_level = 2)
     end
end

# Convergence test with LLF fluxes
open(joinpath(@__DIR__,"log_llf_n2.txt"), "w") do io
     redirect_stdout(io) do
         convergence_test(joinpath(@__DIR__,"elixir_convergence_llf.jl"), 4, polydeg=2, initial_refinement_level = 4)
     end
end
open(joinpath(@__DIR__,"log_llf_n3.txt"), "w") do io
     redirect_stdout(io) do
         convergence_test(joinpath(@__DIR__,"elixir_convergence_llf.jl"), 4, polydeg=3, initial_refinement_level = 3)
     end
end
open(joinpath(@__DIR__,"log_llf_n4.txt"), "w") do io
     redirect_stdout(io) do
         convergence_test(joinpath(@__DIR__,"elixir_convergence_llf.jl"), 4, polydeg=4, initial_refinement_level = 3)
     end
end
open(joinpath(@__DIR__,"log_llf_n5.txt"), "w") do io
     redirect_stdout(io) do
         convergence_test(joinpath(@__DIR__,"elixir_convergence_llf.jl"), 4, polydeg=5, initial_refinement_level = 2)
     end
end


# Convergence test with ES fluxes
open(joinpath(@__DIR__,"log_es_hmat_n2.txt"), "w") do io
    redirect_stdout(io) do
        convergence_test(joinpath(@__DIR__,"elixir_convergence_es_Hmat.jl"), 4, polydeg=2, initial_refinement_level = 4)
    end
end
open(joinpath(@__DIR__,"log_es_hmat_n3.txt"), "w") do io
    redirect_stdout(io) do
        convergence_test(joinpath(@__DIR__,"elixir_convergence_es_Hmat.jl"), 4, polydeg=3, initial_refinement_level = 3)
    end
end
open(joinpath(@__DIR__,"log_es_hmat_n4.txt"), "w") do io
    redirect_stdout(io) do
        convergence_test(joinpath(@__DIR__,"elixir_convergence_es_Hmat.jl"), 4, polydeg=4, initial_refinement_level = 3)
    end
end
open(joinpath(@__DIR__,"log_es_hmat_n5.txt"), "w") do io
    redirect_stdout(io) do
        convergence_test(joinpath(@__DIR__,"elixir_convergence_es_Hmat.jl"), 4, polydeg=5, initial_refinement_level = 2)
    end
end
