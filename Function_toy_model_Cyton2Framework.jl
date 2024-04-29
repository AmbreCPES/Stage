

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
    s::Int = 0
    n_steps::Int = 2000
    activation_rate::Float64 = 0.007
    adata::Vector{Symbol} = [:type, :lineage]
    deaths::DataFrame = DataFrame(id = [], type = [], lineage = [])
    nbr_state::Integer = 6
    death_file::String =  ""
    matrix::Matrix = zeros(nbr_state, nbr_state) 
    pool::Int = 100
    ttnd_parameters::Vector{Tuple{Float64, Float64}} = []
    ttnt_parameters::Vector{Vector{Float64}} = []
    ttd_parameters::Vector{Tuple{Float64, Float64}} =  []

    var_ttnd::Vector{Float64} = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005]
    var_ttnt::Vector{Float64} = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005]
    var_ttd::Vector{Float64} = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005]
end

"""
    HematopoeiticCell

    type::Int : A number corresponding to type of the cell 
    ttd::Int : Time to death in number of steps
    ttnd::Int : Time to next division in number of steps
    ttnt::Int : time to next transition in number of steps

    state::String : "fully_differentiated", "differentiation", "stem"
    time_step_next_event::Int :
    gen::Int : Number of divisions the cell underwent since it left HSC
    lineage::Array{Int, 1} : Array of integer of the form [id, step, type]

"""

@agent HematopoeiticCell NoSpaceAgent begin
    type::Int 
    ttd::Real
    ttnd::Real
    ttnt::Real


    state::String 
    time_next_event::Real #on va comparer le floor de ce nombre de la cellule au time step du model pour savoir si cett cellule doit subir un évènement à ce time step (permet de pas considérer toute les cellules à chaque model step)  
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

function var_distrib(parameters::Float64, n::Int)

        var = rand(Normal(0, parameters), n)
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

function draw_ttd(distribution, parameters::Tuple{Float64, Float64}, n::Int)
    if distribution == LogNormal
       mu = μ_for_mean(parameters[1], parameters[2])
    else 
        mu = parameters[1]
    end
    ttd = rand(distribution(mu, parameters[2]), n)[1]

    #here we are going to draw from a distribution we first suppose lognormal
    # ttnd is set at each division, we can assume it is type dependant
    return ttd
end

#=
HematopoeiticCell properties

type::Int, ttd::Int, ttnd::Int, dd::Int, state::String, time_step_next_event::Int, global_age::Int, gen::Int, lineage::Array{Int, 1} 

=#


#à chaque tour on ne retire pas les paramètres des cellules, quand on réeffectue un tirage, on recompare les 3 horloges, c'est le "next_event"
#dans la life step, si le temps de next division est le plus petit 
#alors, on utilise division2! 
#if ttnd - cell.age < 1 
function division2!(model::ABM, cell::HematopoeiticCell, event_time::Real)
    # un évènement de division se produit 
    
    ttd = var_distrib(model.var_ttd[cell.type], 2) .+ cell.ttd #Variation time to death, we should draw 2 one for each daughter cell
    ttnt = var_distrib(model.var_ttnt[cell.type], 2) .+ cell.ttnt #Variation time to next transition, we should draw 2 one for each daughter cell
    ttnd = draw_ttnd(LogNormal, model.ttnd_parameters[cell.type], 2) .+ model.s .+ event_time  #next division time ( on veut le temps absolu  ou se produit la prochaine division, la distribution determine le temps entre 2 divisions)
    lineage = vcat(cell.lineage, [cell.id, model.s, cell.type])
    gen = cell.gen + 1

    id_daughter_1 = nextid(model)
    id_daughter = [id_daughter_1, id_daughter_1 + 1]
    
    # on regarde quelle sera le prochain évènement des cellules filles
    next_event_1 = findmin([ttd[1] - event_time, ttnd[1], ttnt[1] - event_time]) # la fonction findmin renvoie un tuple (valeur du minimum, indice du minimum)
    next_event_2 = findmin([ttd[2] - event_time, ttnd[2], ttnt[2] - event_time])

    cell_state = [next_event_1[2], next_event_2[2]]

    for event in eachindex(cell_state)
        #on regarde la valeur de next event en fonction on va choisir 
        if cell_state[event] == 1
            add_agent!(
            HematopoeiticCell(id_daughter[event], cell.type, ttd[event], ttnd[event], ttnt[event],"Death" , ttd[event], gen, lineage),
            model
        )

        elseif cell_state[event] == 2
            add_agent!(
            HematopoeiticCell(id_daughter[event], cell.type, ttd[event], ttnd[event], ttnt[event],"Division" ,ttnd[event] , gen, lineage),
            model
        )

        else
            add_agent!(
            HematopoeiticCell(id_daughter[event], cell.type, ttd[event], ttnd[event], ttnt[event], "Transition", ttnt[event], gen, lineage),
            model
        )
        end

    end

    remove_agent!(cell, model)
        
    #on regarde si une des cellules filles a un évènement qui se produit avant la fin de cette time step (se produit dans un temps égal à 1 - event_time)
    #print("next event time =  ", next_event_1[1], "  1 - event_time = ", 1 - event_time, "\n")
    if next_event_1[1] < model.s + 1
        #print("sub life step division2! cell1 \n" )
        life_step!(model[id_daughter[1]], model)
    end

    #print("next event time =  ", next_event_2[1], "  1 - event_time = ", 1 - event_time, "\n")
    if next_event_2[1] < model.s + 1
        #print("sub life step division2! cell2 \n")
        life_step!(model[id_daughter[2]], model)
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
function transition_func(matrix::Matrix, type::Int)
    
    matrix[type,type] = 0

    if sum(matrix[type, : ]) == 0
        new_type = type
    else
        probability_vector = matrix[type, : ]./sum(matrix[type, : ])
        new_type = rand(Categorical(probability_vector), 1)[1]
    end

    return new_type
    
end

#On va écrire un fonction qui permet de tirer selon une distribution binomial un temps jusqu'à la prochaine transition
#Si on cellule dans le type i, on prend la ième ligne de la matrice, le ième coefficient = 1-p (probabilité d'échec) et p probabilité d'une transition vers un autre état. 
#On s'intéresse à la distribution de probabilité d'un succès se produisant au temps t après n échec. 

function transition_time_distribution(matrix_line::Vector{Float64}, cell_type::Int)

    #This is actually a geometric distribution!
    if matrix_line[cell_type] == 1.0
        return [0]
    else
        p = sum(matrix_line) - matrix_line[cell_type]
        print("p = ", p, "\n")
        distrib = [p*((1-p)^n) for n in range(start=0, stop=100, step=0.001)]
        distrib = distrib./sum(distrib)

        cd = [sum(distrib[1:i]) for i in 1:length(distrib)]

        return cd[1:findfirst(x -> x > 0.999, cd)]
    end
end


function draw_ttnt(transition_distribution::Vector{Float64}, step::Int)
    
    random_numbers = rand(1)[1]

    index = findfirst(x -> x > random_numbers, transition_distribution)

    if typeof(index) == Nothing
        print("\n not possible to draw ttnt: transition distribution = ", transition_distribution)
        index = 10000000000 #just a very big number  => transition never happens
    end

    return index/step
end



"""
    transition!(cell::HematopoeiticCell: cell which is going to transition, 
    model : model of the cell, 
    matrix : transition matrix
    )

The cell changes type according to the transition_func(cell, matrix, model.nbr_state) function.
The transition is recorded to the lineage of the cell (cell.lineage, cell.id, model.s, type are added to lineage)

"""

#se produit quand le ttnt - global age est inférieur à 1 à cette time step
function transition2!(model::ABM, cell::HematopoeiticCell, event_time::Real, step::Int)
    
    type = transition_func(model.matrix, cell.type)

    lineage = vcat(cell.lineage, [cell.id, model.s, cell.type])

    cell.lineage = lineage
    cell.type = type

    #cell.ttnd doesn't change it's the time which change so we compare the other time to this one minus event_time
    #but if this event is the next to happend it's time of occurence doesn't change
    #for the other one we had event time because the event of drawing a  new number does not occur at model.s but at model.s + event_time !!! (and what we draw is the time between now and the next event)
    cell.ttnt = draw_ttnt(model.ttnt_parameters[cell.type], step) .+ model.s .+ event_time #ttnt c'est le time step du modele où doit se produire la prochaine transition.

    cell.ttd = draw_ttd(LogNormal, model.ttd_parameters[cell.type], 1) .+ model.s .+ event_time
    next_event = findmin([cell.ttd, cell.ttnd - event_time, cell.ttnt]) 

    if next_event[2] == 1
        cell.state = "Death"
        cell.time_next_event = cell.ttd

    elseif next_event[2] == 2
        cell.state = "Division"
        cell.time_next_event = cell.ttnd

    else
        cell.state = "Transition"
        cell.time_next_event = cell.ttnt
    end

    #on regarde si le prochain éveneemnt se produit au cours de ce time step
    print("cell id = ", cell.id, " next event time =  ", cell.time_next_event)
    if cell.time_next_event < model.s + 1
        life_step!(cell, model)
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

    CSV.write(model.death_file, model.deaths, append = true) #the data of cell which died during this time step is written in a file
    empty!(model.deaths) #dataframe which stores death during one time step is emptied
    model.s += 1 #increase of one of the time step

    return nothing
end

"""
    life_step!(
        cell::HematopoeiticCell : type of agent defined with Agents.jl, 
        model : model defined with Agents.jl 
        )

    Define how an agent evolves to the next time step. 
    life_step is applied to cells which will undergo one of three event during this time step: Division, Transition or Death

"""
function life_step!(cell::HematopoeiticCell, model::ABM)

    if cell.state == "Division"
        print("cell id = ", cell.id, "event time = ", cell.ttnd - model.s,"\n")
        division2!(model, cell, cell.ttnd - model.s) # on fait cell.ttnd - cell.global_age - 1 car on ajoute 1 à global age juste avant

    elseif cell.state == "Transition"
        print("event time = ", cell.ttnt - model.s,"\n")
        transition2!(model, cell, cell.ttnt - model.s, 1000)
    
    elseif cell.state == "Death"
        print("cell id = ", cell.id, "event time = ", cell.ttd - model.s,"\n")
        death!(model, cell, model.deaths)

    end
    return cell.id
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
        n : a function, the model runs until it returns true 
        )

    This methods is a modified version of the run! method defined in Agents.jl, we define this function because we don't need to collect information at each time step.
    We collect information about each agent alive at the last time step.


"""

function custom_run!(model::ABM, agent_step!::Function, model_step!::Function, n::Function)

    step!(model, agent_step!, model_step!, n)
    print("Number of steps: ", model.s, "\n")
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
# Function in building 
##########################################################################################################################################



function ms(model::ABM)

    ids = collect(allids(model))
        # filter all ids whose agents have `w` less than some amount
    filter!(id -> floor(model[id].time_next_event) == model.s, ids)
    return ids
end

function n(model, s)
    while s < model.n_steps

        if length(allids(model)) == 0
            return true
        else
            return false
        end

    end
    return true

end