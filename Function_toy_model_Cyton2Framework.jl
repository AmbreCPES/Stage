

using Agents, DataFrames, CSV, Distributions, StatsBase

#=
I found the Cyton2 model: A model of lymphocyte population and proliferation dynamics, Cheon et al., Front. Bioinform., 2021  (Code published on github: https://github.com/hodgkinlab/cyton2-paper/tree/master) 
=#
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
    adata::Vector{Symbol} = [:type, :lineage]
    deaths::DataFrame = DataFrame(id = [], type = [], lineage = [])
    nbr_state::Integer = 6
    death_file::String =  ""
    matrix::Matrix = zeros(nbr_state, nbr_state) 
    pool::Int = 100
    ttnd_parameters::Vector{Tuple{Float64, Float64}} = []
    ttd_parameters::Vector{Tuple{Float64, Float64}} =  []
end

"""
    HematopoeiticCell

    type::Int : A number corresponding to type of the cell 
    ttd::Int : Time to death in number of steps
    ttnd::Int : Time to next division in number of steps
    ttnt::Int : time to next transition in number of steps

    stage::String : "fully_differentiated", "differentiation", "stem"
    age::Int : age of the cells since the last division
    global_age::Int : number of steps since the cell left the HSC stage
    gen::Int : Number of divisions the cell underwent since it left HSC
    lineage::Array{Int, 1} : Array of integer of the form [id, step, type]

"""

@agent HematopoeiticCell NoSpaceAgent begin
    type::Int 
    ttd::Real
    ttnd::Real
    ttnt::Real


    stage::String 
    time_step_next_event::Int #on va comparer ce nombre de la cellule au time step du model pour savoir si cett cellule doit subir un évènement à ce time step (permet de pas considérer toute les cellules à chaque model step)
    age::Int
    global_age::Int  
    gen::Int 
    lineage::Array{Int, 1} 
    
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
        activation_rate : probability for a cell to become activated at one time step, 
        pool::Vector{Vector} : a collection of cell proprierty 
        )

Add a random number of agent based on the Binomial distribution of parameters (length(pool), activation_rate)

"""
function cells_creation!(model::ABM, pool::Vector{Vector}, activation_rate::Float64)

    activation = rand(Binomial(length(pool), activation_rate), 1)

    for i in 1:activation[1] 

        add_agent!(HematopoeiticCell(nextid(model), pool[i][1], pool[i][2], pool[i][3], pool[i][4], 0, 0, 0, []), model)

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

function death!(model::ABM, cell::HematopoeiticCell, deaths::DataFrame)
    push!(deaths, [cell.id, cell.type, cell.lineage])  
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

function var_ttd(distribution, parameters::Tuple{Float64, Float64}, n::Int)

        var = rand(distribution(parameters[1], parameters[2]), n)
        #here we are probably going to draw from a gaussian centred on 0 .
        #ttd varies at each division 
    return var
end

function draw_ttnd(distribution, parameters::Tuple{Float64, Float64}, n::Int)
    if distribution == LogNormal
       mu = μ_for_mean(parameters[1], parameters[2])
    else 
        mu = parameters[1]
    end
    ttnd = rand(distribution(mu, parameters[2]), n)

    #here we are going to draw from a distribution we first suppose lognormal
    # ttnd is set at each division, we can assume it is type dependant
    return ttnd
end

function number_division(cell::HematopoeiticCell, distribution, )

end
#=
HematopoeiticCell properties

type::Int, ttd::Int, ttnd::Int, dd::Int, stage::String, age::Int, global_age::Int, gen::Int, lineage::Array{Int, 1} 

=#
function division!(model::ABM, cell::HematopoeiticCell, nbr_div::Int) #if the cell divide her generation increase of 1
    mother_cell = cell
    id = nextid(model)
    ids = []
    for i in 1:nbr_div
        push!(ids, range(start = id, stop = (id + 2^(i) - 1), step = 1))
        id += 2^(i)
    end

    for div in eachindex(ids)
        id_mother_cell = 1
        for new_cell in range(start = 1, stop =length(ids[div]), step = 2)
            if div > 1
                mother_cell = model[ids[div - 1][id_mother_cell]]
            end
            v_ttd = var_ttd(Normal, model.ttd_parameters[mother_cell.type], 2) #Variation time to death, we should draw 2 one for each daughter cell
            ttnd = draw_ttnd(LogNormal, model.ttnd_parameters[mother_cell.type], 2) #next division time
            lineage = vcat(mother_cell.lineage, [mother_cell.id, life.s, mother_cell.type])
            gen = mother_cell.gen + 1

        add_agent!(
            HematopoeiticCell(ids[div][new_cell], mother_cell.type, mother_cell.ttd + v_ttd[1], ttnd[1], mother_cell.dd, mother_cell.stage, 0, mother_cell.global_age, gen, lineage),
            model
            )

        add_agent!(
            HematopoeiticCell(ids[div][new_cell + 1], mother_cell.type, mother_cell.ttd + v_ttd[2], ttnd[2], mother_cell.dd, mother_cell.stage, 0, mother_cell.global_age, gen, lineage),
            model
            )

        remove_agent!(mother_cell, model)
        id_mother_cell += 1
        end
    
    end
    return nothing
end

#à chaque tour on ne retire pas les paramètres des cellules, quand on réeffectue un tirage, on recompare les 3 horloges, c'est le "next_event"
#dans la life step, si le temps de next division est le plus petit 
#alors, on utilise division2! 
#if ttnd - cell.age < 1 
function division2!(model::ABM, cell::HematopoeiticCell, event_time::Real)
     
        # un évènement de division se produit 
    v_ttd = var_ttd(Normal, model.ttd_parameters[cell.type], 2) #Variation time to death, we should draw 2 one for each daughter cell
    v_ttnt = var_ttnt(Normal, model.ttnt_parameters[cell.type], 2) #Variation time to next transition, we should draw 2 one for each daughter cell
    ttnd = draw_ttnd(LogNormal, model.ttnd_parameters[cell.type], 2) #next division time
    lineage = vcat(cell.lineage, [cell.id, model.s, cell.type])
    gen = cell.gen + 1

    id_daughter_1 = nextid(model)
    id_daughter_2 = id_daughter_1 + 1
    
    # on regarde quelle sera le prochain évènement des cellules filles
    next_event_1 = findmin([cell.ttd + v_ttd[1] - cell.global_age, ttnd[1], cell.ttnt + v_ttnt[1] - cell.global_age])[2] # la fonction findmin renvoie un tuple (valeur du minimum, indice du minimum)
    next_event_2 = findmin([cell.ttd + v_ttd[2] - cell.global_age, ttnd[2], cell.ttnt + v_ttnt[2] - cell.global_age])[2] 

    cell_stage = [next_event_1, next_event_2]
    for event in eachindex(cell_stage)
        #on regarde la valeur de next event en fonction on va choisir 
        if cell_stage[event] == 1
            cell_stage[event] = "Death"

        elseif cell_stage[event] == 2
            cell_stage[event] = "Division"

        else
            cell_stage[event] = "Transition"
        end
        

    end

    add_agent!(
            HematopoeiticCell(id_daughter_1, cell.type, cell.ttd + v_ttd[1], ttnd[1], cell.ttnt + v_ttnt, cell_stage[1], 0, cell.global_age, gen, lineage),
            model
        )


    add_agent!(
            HematopoeiticCell(id_daughter_2, cell.type, cell.ttd + v_ttd[2], ttnd[2], cell.ttnt + v_ttnt, cell_stage[2], 0, cell.global_age, gen, lineage),
            model
        )

    remove_agent!(cell, model)
        
    #on regarde si une des cellules filles a un évènement qui se produit avant la fin de cette time step (se produit dans un temps égal à 1 - event_time)
    if next_event_1[1] < (1 - event_time)
        life_step!(model[id_daughter_1], model)
    end

    if next_event_2[1] < (1 - event_time)
        life_step!(model[id_daughter_2], model)
    end

    return nothing
end
#end


"""
    transition_func(
        cell::HematopoeiticCell : the cell which is going to transition, 
        matrix::Matrix{Float64} : the transition matrix, 
        nbr_state::Int : the number of possible states for a cell
        )

Draws a random number according to the transition matrix and the type of the cell

"""
function transition_func(cell::HematopoeiticCell, matrix::Matrix{Float64}, nbr_state::Int)

    type = sample(1:nbr_state, ProbabilityWeights(matrix[cell.type, : ]))

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
function transition!(cell::HematopoeiticCell, model::ABM, matrix::Matrix{Float64}, nbr_state::Int)
    type = transition_func(cell, matrix, nbr_state)
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
function model_step!(model::ABM)
    CSV.write(model.death_file, model.deaths, append = true)
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
function life_step!(cell::HematopoeiticCell, model::ABM)
    #cell.age += 1 peut etre qu'on va plutot l'ajouter si cell ne se divise pas soit si celle transitione ou si pas d'évenement
    cell.global_age += 1

    if cell.stage == "Dividing"
    
    elseif cell.stage == "Transition"
    
    elseif cell.stage == "Death"
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
function custom_collect_agent_data!(model::ABM, properties::Vector)
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

function custom_run!(model::ABM, agent_step!::Function, model_step!::Function, n_step::Int)

    step!(model, agent_step!, model_step!, n_step)

    return custom_collect_agent_data!(model, model.adata)
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
function get_adj_list_all_cells(data::DataFrame, last_id::Int, death_file::String)
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
            print(i)
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



function get_adj_list_all_cells2(data, last_id, death_file)
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
    
        

        for i in range(start = L_cell - 2, stop = 4, step = -3)
            print("  ", i)
            if lineage_data[i] != lineage_data[i - 3] 
                
                if (lineage_data[i] in adj_list[lineage_data[i - 3]]) == false
                    
                    push!(adj_list[lineage_data[i - 3]], lineage_data[i])
                    division_type[lineage_data[i - 3]] = lineage_data[i - 1]
                else
                    break
                end
            end
        end
        print("\n")

        if lineage_data[L_cell - 2] != cell_data[cell, :id]
            push!(adj_list[lineage_data[L_cell - 2]], cell_data[cell, :id])
            division_type[lineage_data[L_cell - 2]] = lineage_data[L_cell]
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

function collection_from_lineage(data::DataFrame, last_id::Int, death_file::String)
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


function initialize_model(t0_cells::Vector{Vector{Int}}, ModelParameters) # We expect t0_cells to be an a vector of arrays , each array is a cell with given parameters (type, age, nbr_divisions)
    space = nothing 
    life = ABM(HematopoeiticCell, space; properties = ModelParameters())
    first_cells = HematopoeiticCell.(
                                    1:length(t0_cells), 
                                    map(x -> x[1], t0_cells), 
                                    0, 
                                    0, 
                                    [round(shifted_exponential_sampling(life.division_param[1], life.division_param[2])) for _ in eachindex(t0_cells)] , 
                                    [[i, t0_cells[i][1], 0] for i in eachindex(t0_cells)], 
                                    [round(shifted_exponential_sampling(life.death_param[1], life.death_param[2])) for _ in 1:length(t0_cells)]
                                    )
    
    for cell in  first_cells
        add_agent!(cell, life)
    end

    CSV.write(life.death_file, life.deaths)

    return life
end



##########################################################################################################################################
# Function in building 
##########################################################################################################################################

function death_time_distribution(cell::HematopoeiticCell, model::ABM, mu::Float64, sigma::Float64)

    mu = cell.division_

    return death_time
end