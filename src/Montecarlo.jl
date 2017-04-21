module Montecarlo

using Compat 

export mc_ising_pm1, mc_ising_01, mc_potts 
include("io.jl")
include("ising01_metropolis.jl")
include("isingpm1_metropolis.jl")
include("potts_metropolis.jl")
end
