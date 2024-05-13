

using Agents, DataFrames, CSV, Distributions, StatsBase

#=
I found the Cyton2 model: A model of lymphocyte population and proliferation dynamics, Cheon et al., Front. Bioinform., 2021  (Code published on github: https://github.com/hodgkinlab/cyton2-paper/tree/master) 
=#
#=
This code contains function for a multivariate branching process.

We will first use this process to model Lymphocyte differentiation
=#

##########################################################################################################################################
# Model and agents parameters
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
    n_steps::Int = 40
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
    lineage::Array{Float64, 1} 
    
end


#On va écrire un fonction qui permet de tirer selon une distribution binomial un temps jusqu'à la prochaine transition
#Si on cellule dans le type i, on prend la ième ligne de la matrice, le ième coefficient = 1-p (probabilité d'échec) et p probabilité d'une transition vers un autre état. 
#On s'intéresse à la distribution de probabilité d'un succès se produisant au temps t après n échec. 
"""
    transition_time_distribution(
        matrix_line::Vector{Float64} : line in the transition matrix corresponding to the type of the cell, 
        cell_type::Int : type of the cell.

        stop::Int = 100 : the maximum number of step before transition, 
        step::Float64 = 0.001 : the time between 2 points in the cumulative distribution
        )

Setting of the transition time discrete cumulative distribution given with the transition matrix.
Probability distribution corresponds to the probability a transition occurs for the first time. probability to transition is equal to 1 - probability to stay in the same state.

"""
function transition_time_distribution(matrix_line::Vector{Float64}, cell_type::Int; stop::Int = 100, step::Float64 = 0.001)

    if matrix_line[cell_type] == 1.0
        return [0]
    else
        p = sum(matrix_line) - matrix_line[cell_type]
        distrib = [p*((1-p)^n) for n in range(start=0, stop=stop, step=step)]
        distrib = distrib./sum(distrib)

        cd = [sum(distrib[1:i]) for i in 1:length(distrib)]

        return cd[1:findfirst(x -> x > 0.999, cd)]
    end
end


##########################################################################################################################################
# Functions to evolve agents at each step
##########################################################################################################################################

#We define function to draw time to death, time to transition and time to division but also to perturb these time

"""
    var_distrib(
        parameters::Float64 : standart deviation of the gaussian given the cell type, 
        n::Int : number of number generated
    )

This methods allow to add a Gaussian perturbation to a Time (ttd, ttnd or ttnt).
"""
function var_distrib(parameters::Float64, n::Int)

    var = rand(Normal(0, parameters), n)
    #here we are probably going to draw from a gaussian centred on 0 .
    #ttd varies at each division 
return var
end


"""
    draw_ttnd(
        distribution: a continuous distribution defined in Distributions.jl (often LogNormal), 
        parameters::Tuple{Float64, Float64} : (mu, sigma) the parameters of the distribution (parameters are cell type dependant), 
        n::Int : number of ttnd we wish to draw from this distribution
        )

Function to draw a new time to next division. We draw the time between to division given a cell type.
New ttnd are drawn at each division.

"""
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

"""
    draw_ttd(
        distribution: a continuous distribution defined in Distributions.jl (often LogNormal), 
        parameters::Tuple{Float64, Float64} : (mu, sigma) the parameters of the distribution (parameters are cell type dependant), 
        n::Int : number of ttnd we wish to draw from this distribution
        )

Function to draw a new time to death. We draw the time to death for a cell type.
A new ttd is drawn at each transition.

"""
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
"""
    draw_ttnt(
        transition_distribution::Vector{Float64} : discrete cumulative distribution of the transition time for a cell type, 
        step::Int : the time between 2 points in the cumulative distribution
        )

Function to draw a new time to next transition. We draw the time to transition for a cell type.
A new ttnt is drawn at each transition.

"""
function draw_ttnt(transition_distribution::Vector{Float64}, step::Int)
    
    random_numbers = rand(1)[1]

    index = findfirst(x -> x > random_numbers, transition_distribution)

    if typeof(index) == Nothing
        index = 10000000000 #just a very big number  => transition never happens
    end

    return index/step
end

function draw_ttnt2(matrix_line::Vector{Float64}, cell_type::Int)

    if matrix_line[cell_type] == 1.0
        return [10000]
    else
        p = 1 - matrix_line[cell_type] #probability of success
        #Initialization exponential distriution 
        theta = −(log(1 − p))
        x =rand(1)[1]
        
       ttnt = -log(1 - x)/theta 

        return [ttnt, x]
    end
end

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


#We define function for cell event: death, division, transition

"""
    death!(
        cell::HematopoeiticCell : the cell which will be removed from the model,
        model:: ABM : model from which the cell is removed
        deaths::DataFrame : DataFrame to which the parameters of the dead cells will be added
    )

Remove cell from the model and record these information: [cell.id, cell.type, cell.lineage] in a DataFrame

"""

function death!(model::ABM, cell::HematopoeiticCell, deaths::DataFrame)
    cell.lineage = vcat(cell.lineage, [cell.id, round(cell.ttd, digits = 3), cell.type])
    push!(deaths, [cell.id, cell.type, cell.lineage])  
    remove_agent!(cell, model)

	return nothing
end


#=
HematopoeiticCell properties

type::Int, ttd::Int, ttnd::Int, dd::Int, state::String, time_step_next_event::Int, global_age::Int, gen::Int, lineage::Array{Int, 1} 

=#
"""
    division2!(
        model::ABM : model in which the cell divides, 
        cell::HematopoeiticCell : cell which divides, 
        event_time::Real : Float64 in [0;1], the time in the time step when the event occurs
    )

2 daughter cell are created. They inherit their ttd, ttnt from the mother cells (+- gaussian perturbation) and are attributed a new time to next division which differs
between the daughter from gaussian perturbation.

We determine the next event of the cells: "Death", "Division", "Transition" and their time.
If the next event is due to happen before the next time step (cell.time_next_event < model.s + 1) then the cell undergoes a life_step!

"""
function division2!(model::ABM, cell::HematopoeiticCell, event_time::Real)
    # un évènement de division se produit 
    
    ttd = var_distrib(model.var_ttd[cell.type], 2) .+ cell.ttd #Variation time to death, we should draw 2 one for each daughter cell
    ttnt = var_distrib(model.var_ttnt[cell.type], 2) .+ cell.ttnt #Variation time to next transition, we should draw 2 one for each daughter cell
    ttnd = draw_ttnd(LogNormal, model.ttnd_parameters[cell.type], 2) .+ model.s .+ event_time  #next division time ( on veut le temps absolu  ou se produit la prochaine division, la distribution determine le temps entre 2 divisions)
    lineage = [vcat(cell.lineage, [cell.id, round(cell.ttnd, digits = 3), cell.type]), [cell.id, round(cell.ttnd, digits = 3), cell.type]]
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
            if ttd[event] < cell.ttnd
                ttd[event] = cell.ttnd
            end
            add_agent!(
            HematopoeiticCell(id_daughter[event], cell.type, ttd[event], ttnd[event], ttnt[event],"Death" , ttd[event], gen, lineage[event]),
            model
        )

        elseif cell_state[event] == 2
            if ttnd[event] < cell.ttnd
                ttnd[event] = cell.ttnd
            end
            add_agent!(
            HematopoeiticCell(id_daughter[event], cell.type, ttd[event], ttnd[event], ttnt[event],"Division" ,ttnd[event] , gen, lineage[event]),
            model
        )

        else
            if ttnt[event] < cell.ttnd
                ttnt[event] = cell.ttnd
            end
            add_agent!(
            HematopoeiticCell(id_daughter[event], cell.type, ttd[event], ttnd[event], ttnt[event], "Transition", ttnt[event], gen, lineage[event]),
            model
        )
        end

    end

    remove_agent!(cell, model)
        
    if next_event_1[1] < model.s + 1
        life_step!(model[id_daughter[1]], model)
    end

    if next_event_2[1] < model.s + 1
        life_step!(model[id_daughter[2]], model)
    end

    return nothing
end
#end


"""
    transition_func(
        matrix::Matrix{Float64} : the transition matrix of the model, 
        type::Int : type of the cell which is transitionning
        )

Draws a random number according to the transition matrix and the type of the cell

"""
function transition_func(matrix::Matrix, type::Int)
    matrix_line = matrix[type,:]
    matrix_line[type] = 0
    if sum(matrix[type,:]) != 1
        return error("la matrice de transition n'est pas une matrice de probabilité")
    else
        if sum(matrix_line) == 0
            new_type = type
        else
            probability_vector = matrix_line./sum(matrix_line)
            new_type = rand(Categorical(probability_vector), 1)[1]
        end
    end
    return new_type
    
end


"""
    transition!(cell::HematopoeiticCell: cell which is going to transition, 
    model : model of the cell, 
    matrix : transition matrix
    )

The cell changes type according to the transition_func(cell, matrix, model.nbr_state) function.
The transition is recorded to the lineage of the cell (cell.lineage, cell.id, model.s, type are added to lineage)

"""

function transition2!(model::ABM, cell::HematopoeiticCell, event_time::Real)
    
    type = transition_func(model.matrix, cell.type)

    lineage = vcat(cell.lineage, [cell.id, round(cell.ttnt, digits = 3), cell.type])

    cell.lineage = lineage
    cell.type = type

    #cell.ttnd doesn't change it's the time which change so we compare the other time to this one minus event_time
    #but if this event is the next to happend it's time of occurence doesn't change
    #for the other one we had event time because the event of drawing a  new number does not occur at model.s but at model.s + event_time !!! (and what we draw is the time between now and the next event)
    cell.ttnt = draw_ttnt2(model.matrix[cell.type,:], cell.type)[1] + model.s + event_time #ttnt c'est le time step du modele où doit se produire la prochaine transition.
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
    if cell.time_next_event < model.s + 1
        life_step!(cell, model)
    end
    
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
        division2!(model, cell, cell.ttnd - model.s) # on fait cell.ttnd - cell.global_age - 1 car on ajoute 1 à global age juste avant

    elseif cell.state == "Transition"
        transition2!(model, cell, cell.ttnt - model.s)
    
    elseif cell.state == "Death"
        death!(model, cell, model.deaths)

    end
    return cell.id
end





##########################################################################################################################################
# Functions for model evolution
##########################################################################################################################################


"""
    ms(
        model::ABM : an ABM from Agents.jl
    )

A scheduler which determines the cell which will undergo life_step! at each model_step!.
Every cell with an event scheduled during the current time step (floor(cell.time_next_event) == model.s) is scheduled to undergo life_step at the current time step

"""

function ms(model::ABM)

    ids = collect(allids(model))
        # filter all ids whose agents have `w` less than some amount
    filter!(id -> floor(model[id].time_next_event) == model.s, ids)
    return ids
end


"""
    n(
        model::ABM : an ABM from Agents.jl,
        s::Int : the current time step
    )
    
n determines the number of steps the model will run for. While n returns false, the model run.
Here the model runs while there are cells alive. Even if cells are still alive, it will stop at the max time step, model.n_steps
"""
function n(model, s)
    while s < (model.n_steps-1)

        if length(allids(model)) == 0
            return true
        else
            return false
        end

    end
    return true

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
    data = custom_collect_agent_data!(model, model.adata)

    data[!,:lineage] .= [vcat(data[x,:lineage], [data[x, :id], model.s + 1, data[x, :type]]) for x in 1:nrow(data)]

    return data
end

