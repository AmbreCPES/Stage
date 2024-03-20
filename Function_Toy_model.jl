export ModelParameters, HematopoeiticCell, custom_collect_agent_data!, custom_run!, get_adj_list_all_cells, collection_from_lineage

using Agents, DataFrames, CSV, Distributions
 
#=
This code contains function for a multivariate branching process.

We will first use this process to model Lymphocyte differentiation
=#

##########################################################################################################################################
# Defining model parameters and agents type
##########################################################################################################################################
"""
   ModelParameters

    s::Int = 0 : the time step of the model
    activation_rate::Float : rate at which cell are created (temporary as long as niche not implemented)
    deaths::DataFrame 
    nbr_state::Int : for rand()  to generate a matrix and not a vector of vector nbr_state must be a Int64 and not a Int8
    death_file::String : file to register the deaths
    matrix::Matrix : Initialisation of the transition matrix, here at random
    death_param::Tuple{Float, Float} : parameter of the shifted exponential used to sample division time
    division_param::Tuple{Float, Float}
    pool::Int : an integer representing the number of quiescent cells which can be activated

   Example: 
    parameters = ModelParameters(0, 1/700, DataFrame(), 6, "hello", zeros(6,6),(0.1, 0.1),(0.1, 0.1))

"""

Base.@kwdef mutable struct ModelParameters
    s::Int = 1 
    activation_rate::Float64 = 0.007
    deaths::DataFrame = DataFrame(id = [], type = [], lineage = [])
    nbr_state::Integer = 6
    death_file::String =  "C:\\Users\\ambre\\Documents\\ENS\\stage_M1\\code_1\\test_file"
    matrix::Matrix = zeros(nbr_state, nbr_state)
    death_param::Tuple{Float64, Float64} = shifted_exponential_param([0, 80],[0, 0.99])
    division_param::Tuple{Float64, Float64} = shifted_exponential_param((17.8, 57), (0.77, 0.99))  
    pool::Int =100
end

"""
    HematopoeiticCell

    type::Int : A number corresponding to type of the cell 
    age::Int : Age of the cell since she left HSC in number of time step
    divisions::Int : number of division of a cell since they left the stem cell state
    nbr_division::Int : number of divisions until the next division
    lineage::Array{Int, 1} : an array of int of the form [..., id, step, type, ...]
    divsion_time::Int : number of time step since the last division

    Example:
    HematopoeiticCell(1, 1, 0, 0, 15, [], 0) , the first argument is the ID of the cell, a mandatory field of NoSapceAgent

"""

@agent HematopoeiticCell NoSpaceAgent begin
    type::Int
    age::Int 
    divisions::Int 
    nbr_division::Int 
    lineage::Array{Int, 1} 
    divsion_time::Int
end

##########################################################################################################################################
# Defining Sampling Methods
##########################################################################################################################################

"""
    shifted_exponential_param(T,X) 

X an T are tuples of 2 floats: T = (t1, t2) and X = (x1, x2), which corresponds to point coordinates
T is the absciss and X the ordonate.

This function calculate the parameter of the  shifted exponential cumulative distribution given 2 points.

return the parameter lambda and L as a tuple.

Example:
    shifted_exponential_param((17.8, 57), (0.77, 0.99))

"""
function shifted_exponential_param(T,X) 
    L = ( T[2] * log(1 - X[1])/(log(1 - X[2])) + T[1] ) * 1/(1 - log(1 - X[1]))
    lambda = log(1 - X[2])/(L - T[2])
    return (lambda, L)
end

"""
    shifted_exponential_sampling(
        Lambda, L : Lambda and L are the parameters of the shifted exponential cumulative distribution
        ) 

This function allow sampling according to shifted exponential of the given parameters.

return a float number 

Example:
    param = shifted_exponential_param((17.8, 57), (0.77, 0.99))
    results = shifted_exponential_sampling(param[1], param[2])

"""

function shifted_exponential_sampling(lambda, L)
     alea = rand(Float16)
     return  -log(1 - alea)/lambda + L
 end


##########################################################################################################################################
# Defining functions to evolve model and agents each step
##########################################################################################################################################

###########
# TO TEST
###########
"""
    cells_creation!(
        model : the model where the cells will be added,
        pool::Integer : number of quiescent cell which can become activated,
        )

Add a random number of agent based on the Binomial distribution of parameters (pool, activation_rate)

"""
function cells_creation!(model, pool)
    activation_rate = model.activation_rate
    activation = rand(Binomial(pool, activation_rate), 1)

    for _ in 1:activation[1] 
        
        add_agent!(HematopoeiticCell(nextid(model), 1, 0, 0, 0, [], round(shifted_exponential_sampling(model.division_param[1], model.division_param[2]))), model)
    end

    return nothing
end


"""
    death!(
        cell::HematopoeiticCell : the cell which will be removed from the model,
        model: the model from which the cell is drawn
    )

Remove cell from the model and record these information: [cell.id, cell.type, cell.lineage] in a DataFrame

"""
function death!(cell::HematopoeiticCell, model)
    push!(model.deaths, [cell.id, cell.type, cell.lineage])  
    remove_agent!(cell, model)

	return nothing
end

"""
    division!(
        cell::HematopoeiticCell : a cell which is dividing,
        model : the model where the cell is
        )

Create two cells with parameters of the dividing cell and remove the dividing cell from the model.
We add [cell.id, model.s, cell.type] and the lineage of the mother cell and the sum is given to daughter cells.
The dividing time of each daughter cell is drawn from shifted exponential distribution.

"""

function division!(cell::HematopoeiticCell, model)
    divisions = cell.divisions + 1
    lineage = vcat(cell.lineage, [cell.id, model.s, cell.type])
    add_agent!(HematopoeiticCell(nextid(model), cell.type, 0, divisions, cell.nbr_division, lineage, round(shifted_exponential_sampling(model.param[1], model.param[2]))), model)
    add_agent!(HematopoeiticCell(nextid(model), cell.type, 0, divisions, cell.nbr_division, copy(lineage), round(shifted_exponential_sampling(model.param[1], model.param[2]))), model)
    remove_agent!(cell, model)

    return nothing
end

"""
    transition_func(
        cell::HematopoeiticCell : the cell which is going to transition, 
        matrix::Matrix{Float64} : the transition matrix, 
        nbr_state::Int : the number of possible states for a cell
        )

Draws a random number according to the transition matrix and the type of the cell

"""
function transition_func(cell::HematopoeiticCell, matrix::Matrix{Float64}, nbr_state::Int)
    type = sample(1:nbr_state, ProbabilityWeights(matrix[ : ,cell.type]))

    return type
end

"""
    transition!(cell::HematopoeiticCell: cell which is going to transition, 
    model : model of the cell, 
    matrix : transition matrix
    )

The cell changes type according to the transition_func(cell, matrix, model.nbr_state) function.
The transition is recorded to the lineage of the cell (cell.lineage, cell.id, model.s, type are added to lineage)

"""
function transition!(cell::HematopoeiticCell, model, matrix)
    type = transition_func(cell, matrix, model.nbr_state)
    if type != cell.type
        push!(cell.lineage, cell.id, model.s, type)
        cell.type = type
    end
    
    return nothing
end

"""
    model_step!(
        model : model defined with Agents.jl 
        )

define how the model evolves to the next time step

"""
function model_step!(model)
    CSV.write(model.file, model.deaths, append = true)
    empty!(model.deaths)
    model.matrix = rand(Dirichlet(model.nbr_state, 2), model.nbr_state)
    model.s += 1

    return nothing
end

"""
    life_step!(
        cell::HematopoeiticCell : type of agent defined with Agents.jl, 
        model : model defined with Agents.jl 
        )

    define how an agent evolves to the next time step

"""
function life_step!(cell::HematopoeiticCell, model)
    cell.age += 1
    if cell.divisions != cell.nbr_division
        transition!(cell, model, model.matrix)
        if cell.age == cell.divsion_time
            division!(cell, model)
        end
    else
        death!(cell, model)
    end

    return nothing
end



##########################################################################################################################################
# Defining functions to run the model
##########################################################################################################################################

"""
    custom_collect_agent_data!(
        model : model from which the data are to be collected, 
        properties::Vector : properties to be collected for each agent
        )

Collect properties of all cells after the model is ran
"""
function custom_collect_agent_data!(model, properties::Vector)
    alla = sort!(collect(allagents(model)), by = a -> a.id)
    dd = DataFrame()
    dd[!, :id] = map(a -> a.id, alla)
    for fn in properties
        # Extract data for each property using getfield
        dd[!, fn] = getfield.(alla, fn)
    end

    return dd
end

"""
    custom_run!(
        model: A model as defined in Agents.jl which is evolved by the function,
        agent_step!: a function which determines how each agent will evolve to the next step,
        model_step!: a function which determines how the model will evolve to the next step
        n_step : number of steps to run the model over
        )

    This methods is a modified version of the run! method defined in Agents.jl, we define this function because we don't need to collect information at each time step.
    We collect information about each agent at the model.n_step (model property) time step.


"""

function custom_run!(model, agent_step!, model_step!, n_step)

    step!(model, agent_step!, model_step!, n_step)

    return custom_collect_agent_data!(model, adata)
end


##########################################################################################################################################
# Plotting data as graph
##########################################################################################################################################


###########
# TO TEST
###########

"""
    get_adj_list_all_cells(
        data : DataFrame containing data of cells at the end of the run,
        last_id : last given id,
        death_file : Path of file containing CSV file containing data of the cells which died during the run 
        )

Return adjacent matrix from the data (asymetrix matrix which cooresponds to a directed graph) and type of each cell as a dictionnary (key = cell id and value = type of the cell)

"""

function get_adj_list_all_cells(data, last_id, death_file)
    deaths = DataFrame(CSV.File(death_file))
    cell_data = vcat(data, deaths)
    adj_list = [[] for _ in 1:last_id]
    division_type = Dict{Int, Int}(cell_data[!, :id][i] => cell_data[!, :type][i] for i in 1:nrow(cell_data))
    
    for cell in 1:nrow(cell_data)
        if isa(cell_data[cell, :lineage], String)
            cell_data[cell, :lineage] = parse.(Int, split(chop(cell_data[cell, :lineage], head = 1, tail = 1), ","))
        end
        lineage_data = cell_data[cell, :lineage]
    
        sort!(cell_data, :id)
    
        L_cell = length(lineage_data)
    
        if lineage_data[L_cell - 2] != cell_data[cell, :id]
            push!(adj_list[lineage_data[L_cell - 2]], cell_data[cell, :id])
            division_type[lineage_data[L_cell - 2]] = lineage_data[L_cell]
        end

        for i in range(start = L_cell - 2, stop = 4, step = -3)
            if lineage_data[i] != lineage_data[i - 3]  
                if (lineage_data[i] in adj_list[lineage_data[i - 3]]) == false
                    push!(adj_list[lineage_data[i - 3]], lineage_data[i])
                    division_type[lineage_data[i - 3]] = lineage_data[i - 1]
                end
            end
        end
    end

    return adj_list, division_type
end


###########
# TO TEST
###########
"""
    collection_from_lineage(
        data : DataFrame containing data of cells at the end of the run, 
        deaths : DataFrame containing data of the cells which died during the run , 
        last_id : last given id
        )

Return adjacent matrix from the data (asymetrix matrix which cooresponds to a directed graph) and type of each cell as a dictionnary (key = cell id and value = type of the cell)

"""

function collection_from_lineage(data, last_id, death_file)
    #data is a DataFrame each row corresponds to an alive cell with colunms id, type, lineage
    #last_id is the last ID of all cells which have existed
    #:lineage is in the form of an array recording every transition and division in the form of [id, step, type, id, step, type, ...]
    
    deaths = DataFrame(CSV.File(death_file))
    cell_data = vcat(data, deaths)
    adj_matrix = zeros(Int, last_id, last_id)
    division_type = Dict{Int, Int}(cell_data[!, :id][i] => cell_data[!, :type][i] for i in 1:nrow(cell_data))
    for cell in 1:nrow(cell_data)

        if isa(cell_data[cell, :lineage], String)
            cell_data[cell, :lineage] = parse.(Int, split(chop(cell_data[cell, :lineage], head = 1, tail = 1), ","))
        end
        lineage_data = cell_data[cell, :lineage]
        L = length(lineage_data)

        for i in range(start = 4, stop = L, step = 3)
            if lineage_data[i] != lineage_data[i - 3]
                adj_matrix[lineage_data[i - 3], lineage_data[i]] = 1
                division_type[lineage_data[i - 3]] = lineage_data[i - 1]
            end
        end
    end

    return adj_matrix, division_type
end



##########################################################################################################################################
# Model Initialisation
##########################################################################################################################################


function initialize_model(t0_cells, parameters) # We expect t0_cells to be an a vector of arrays , each array is a cell with given parameters (type, age, nbr_divisions)
    space = nothing 
    life = ABM(HematopoeiticCell, space; properties = parameters)
    first_cells = HematopoeiticCell.(1:length(t0_cells), map(x -> x[1], cell0), 0, 0, [round(shifted_exponential_sampling(life.division_param[1], life.division_param[2])) for _ in 1:length(t0_cells)] , [[] for i in 1:length(t0_cells)], [round(shifted_exponential_sampling(life.death_param[1], life.death_param[2])) for _ in 1:length(t0_cells)])
    
    for cell in  first_cells
        add_agent!(cell, life)
    end

    CSV.write(life.death_file, life.deaths)

    return life
end




