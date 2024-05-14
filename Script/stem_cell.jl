#Stem cell modelling
include("Function_toy_model_Cyton2Framework.jl")


Base.@kwdef mutable struct NicheParameters
    stem_cell_number::Int = 100
    activation_rate::Float64 = 0.007
    K::Int = #Carrying capacity
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
function niche_initialization(nbr_stem_cells::Int, K, division_rate_parameters::Tuple)

end

#2 separates ABM one representing the lineage and one the niche 
function stem_cell_count(model::ABM) #to use in the model step

end

function division_proba()
end

function hematopoeitic_cell_creation()
end

