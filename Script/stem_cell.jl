#Stem cell modelling
include("Function_toy_model_Cyton2Framework.jl")
using Pkg
Pkg.activate("/home/ajaeger/Documents/githubrepos/Stage/Script/MyProject") #activation of the correct environment for the notebook

Base.@kwdef mutable struct NicheParameters
    stem_cell_number::Int = [100]
    activation_parameters::Tuple = Normal(1/110)
    K::Int = 17000 #Carrying capacity
    division_time_parameters::Tuple = Normal(0.0131)
end

@agent StemCell NoSpaceAgent begin
    type::Int 
    ttd::Real
    ttnd::Real
    ttnt::Real

    state::String 
    time_next_event::Real #on va comparer le floor de ce nombre de la cellule au time step du model pour savoir si cett cellule doit subir un évènement à ce time step (permet de pas considérer toute les cellules à chaque model step)  
    gen::Int 
    lineage::Array{Float64, 1} 
end

function division_prob(model::ABM)
    
end
function niche_initialization(nbr_stem_cells::Int, K, division_rate_parameters::Tuple)

end

#2 separates ABM one representing the lineage and one the niche 
function stem_cell_count(model::ABM) #to use in the model step

end

function division_proba()
end

function hematopoeitic_cell_creation()
end

using DifferentialEquations

#Adam mac Lean M3 model 2014 Paper

#rho = 1/110 #growth rate 
#K = 17000   #Carrying capacity
#delta = 0.1 # "death rate" actually its more a rate of "leaving the stem cell state wether through death or differentiation 
#phi = x, rho, K -> (rho + 1)*K/(K+(rho*x))
#function phi(rho, K, x)
 #   return (rho + 1)*K/(K+(rho*x))
#end
#p = [rho, K, delta]
function log_growth!(du,u,p,t) 
    phi = (p[1] + 1)*p[2]/(p[2]+(p[1]*u[1]))
    du[1] = (phi - 1)*p[3]*u[1]
    du[2] = p[3]*u[1]
end

prob = ODEProblem(log_growth!, [1580.0; 0.0], (0.0, 15000.0), p = (0.1, 1700, 0.01))
sol = solve(prob)