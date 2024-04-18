#= 
When cells leave the stem cell state they must be attributed characteristics which will define their fate
=#
include("Function_toy_model_Cyton2Framework.jl")

using SpecialFunctions, Distributions

function sampling_lognorm(mu, sigma, n_cells)

    return rand(LogNormal(mu, sigma), n_cells)
end

#= For each type, we have 3 distributions:
    - time to death
    - time to next division
    - division destiny

We store it in the form of a stuct of properties: distribution's parameters time to death, distribution's parameters time to next division, distribution's parameters division destin
Parameters are in Tuple, the length of the tuple depends on the nature of the distribution
The index of the vector corresponds to its type. 

=#
function μ_for_mean(m, σ)
    return log(m) - σ^2/2
end


function initialize_distributions(typeparameters::TypeParameters, n_states::Int)
    #For now all distributions are supposed to be LogNormal
    type_distribution = [[] for _ in 1:n_states]
    for i in 1:n_states
        type_distribution[i] = [LogNormal(typeparameters.ttd[i][1], typeparameters.ttd[i][2]),
                                LogNormal(typeparameters.ttnd[i][1], typeparameters.ttnd[i][2])]
    end
    return type_distribution
end

#creation of a collection of cell
# As entry a list of tuple of the form collection_t0 = [(type, number_of_cells_of_this_type), ...]
# To each cell we attribut a death time, a number of divisions, a time to next division, we draw them distributions here log normal distribution

#=
HematopoeiticCell properties

type::Int, ttd::Int, ttnd::Int, dd::Int, stage::String, age::Int, global_age::Int, gen::Int, lineage::Array{Int, 1} 

=#

function initialize_collection(type_distributions::Vector{Vector{Any}},transition_distribution::Vector{Vector{Float64}}, collection_t0::Vector{Tuple{Int64, Int64}}, n_tot::Int)
    # we create a vector of the form [[id, type, time_to_death, time_to_next_division, time_to_next_transition], ...]
    cell_collection = [[] for _ in 1:n_tot]
    ni = 0 #the index of the last cell which received properties

    for type in eachindex(collection_t0)
        n_cells = collection_t0[type][2]
        push!.(
            cell_collection[ni + 1 : ni + n_cells], 
            [collection_t0[type][1] for _ in 1:n_cells], 
            rand(type_distributions[type][1], n_cells), 
            rand(type_distributions[type][2], n_cells), 
            [draw_ttnt(transition_distribution[collection_t0[type][1]]) for _ in 1:n_cells]
            )
        ni += n_cells
    end
    return cell_collection
end

function initialize_modelparameters(nbr_state::Int, death_file::String, matrix::Matrix, typeparameters_0::TypeParameters, transition_distribution_0::TransitionDistribution)
    modelparameters_0 = ModelParameters(nbr_state = nbr_state, death_file = death_file, matrix = matrix, ttnd_parameters = typeparameters_0.ttnd, ttnt_parameters = transition_distribution_0.ttnt, ttd_parameters = typeparameters_0.ttd)
    
    return modelparameters_0
end


function initialize_model(collection_t0::Vector{Tuple{Int64, Int64}},transition_distribution::TransitionDistribution, modelparameters_0::ModelParameters, typeparameters_0::TypeParameters, n_tot::Int, ms) # We expect t0_cells to be an a vector of arrays , each array is a cell with given parameters (type, age, nbr_divisions)
    space = nothing

    n_states = length(typeparameters_0.ttd)
    type_distribution_0 = initialize_distributions(typeparameters_0, n_states)
    cell_collection_0 = initialize_collection(type_distribution_0, transition_distribution.ttnt, collection_t0, n_tot)

    life = StandardABM(HematopoeiticCell, space; properties = modelparameters_0, scheduler = ms)

#= type::Int ttd::Real ttnd::Real ttnt::Real   state::String  time_step_next_event::Int  gen::Int  lineage::Array{Int, 1} 
 
=#
    first_cells = HematopoeiticCell.(
                                    1:n_tot,
                                    map(x -> x[1], cell_collection_0),
                                    map(x -> x[2], cell_collection_0), 
                                    map(x -> x[3], cell_collection_0),
                                    map(x -> x[4], cell_collection_0),
                                    ["initialisation" for _ in 1:n_tot],
                                    zeros(Int, n_tot), 
                                    zeros(Int, n_tot),
                                    [[] for _ in 1:n_tot])
    
    for cell in  first_cells
        next_event = findmin([cell.ttd, cell.ttnd, cell.ttnt]) 
        if next_event[2] == 1
            cell.state = "Death"

        elseif next_event[2] == 2
            cell.state = "Division"

        else
            cell.state = "Transition"
        end

        cell.time_step_next_event = floor(next_event[1])
        add_agent!(cell, life)
    end

    CSV.write(life.death_file, life.deaths)

    return life
end


"""
set_transition_distrib(parameters::Vector{Tuple})

    parameters::Vector{Tuple} : vector of the form [(type, mu, sigma), ...]

    this function initialise a Dict of Gaussian distribution given the parameters

"""                     
function set_transition_distrib(parameters::Vector)
    transition_distrib = Dict()
    for i in parameters
        transition_distrib[i[1]] = Normal(i[2], i[3])
    end
    return transition_distrib
end

"""
set_transition_matrix!(
    mat::Matrix{Float64}
    )
    This function modifies the diagonal coefficient so that the sum of each line is equal to one
"""
function set_transition_matrix!(mat::Matrix{Float64})
    size_mat = size(mat)[1]
    for x in range(start = 1, stop = size_mat, step = 1)
        mat[x,x] = 0
        mat[x,x] = 1 - sum(mat[x, : ])
    end

    sum_line = sum_matrix_line(mat, size_mat)

    if sum_line == [1.0 for _ in 1:size_mat]  
        return nothing
    else 
        return sum_line
    end
end

function sum_matrix_line(mat::Matrix, size_mat::Int)
    sum_line = [sum(mat[x, : ])  for x in range(start = 1, stop = size_mat, step = 1)]
    return sum_line
end

"""
update_transition_matrix(
    mat::Matrix{Float64}, 
    distributions::Dict : distribution assciated with a type (of the same length as ttu)
    ttu::Vector{Int} : type to update a vector containing the type to update (line in the matrix), the index in the vector is the line, the vector at this index contains the type to update, vector is empty if there is nothing to update
    )
    
    We draw coefficient for the transition matrix from different distribution for each type which can undergo transition
    
"""

function variation_transition_matrix!(mat::Matrix{Float64}, parameters::Vector, ttu::Vector{Vector{Any}})
    
    for i in 1:length(ttu)
        if length(ttu[i]) != 0
            mat[i, ttu[i]] = mat[i, ttu[i]] + rand(truncated(Normal(0, parameters[i]), -1, 1), length(ttu[i]))
            normalize!(mat[i,:],1)
        end
    end

    return nothing
end