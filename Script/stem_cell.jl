#Stem cell modelling
using Pkg
Pkg.activate("/home/ajaeger/Documents/githubrepos/Stage/Script/MyProject") #activation of the correct environment for the notebook
using Agents, Distributions, DataFrames, DifferentialEquations, Plots, GLM, LinearRegression

##################################################################################################################################
# Model Parameters 
##################################################################################################################################

Base.@kwdef mutable struct LayersNicheParameters
    s::Int = 0
    cell_created::Vector = [] #everytime a cell divide we count that one cell is created 
    cell_differentiated::Vector = [] #everytime a cell differentiate we count + 1 
    stem_cell_number::Vector = [[], [], []] #at index 1 is the number of quiescent stem cell and at index 2 the number of cycling stem cell
    activation_parameters::Vector{Tuple} = [(log(2)/(1/110), 10), (log(2)/(1/110), 10)] #at index 1 is the time of activation of quiescent stem cell and at index 2 the time of differentiation of cycling stem cell
    quiescence_parameters::Vector{Tuple} = [(log(2)/(1/110), 10)] #time of cycling cell becoming quiescent again
    K::Int = 17000 #Carrying capacity
    division_parameters::Vector{Tuple} = [(log(2)/(1/110), 10)] #1/110 is the stem cell proliferation rate
end


@agent ActivableStemCell NoSpaceAgent begin
    type::Int #quiescent stem cell are type 0 and cycling are type 1
    division_time::Float64
    activation_time::Float64
    quiescence_time::Float64
    differentiation_time::Float64
    label::Bool
    next_event_time::Float64
    next_event_type::String #Differentiation, Division, Activation, Quiescence
end


##################################################################################################################################
# Model Functions
##################################################################################################################################

function activation!(model::ABM, cell::ActivableStemCell)

    cell.type = 1
    cell.activation_time = 0.0
    
    quiescence_differentiation!(model, cell, cell.next_event_time)

    if (cell.division_time < model.s) | (cell.quiescence_time < model.s) | (cell.differentiation_time < model.s)
        print("\n Y a un probleme dans l'attribution des temps \n")
        print("On est au temps : ", model.s)
        print("\n l'évenement se deroule au temps : ", cell.next_event_time)
        print("Division time = ", cell.division_time, " Quiescence_time = ", cell.quiescence_time, " differentiation time = ", cell.differentiation_time)
    end

    if cell.next_event_time < model.s + 1
        agent_step!(cell, model)
    end
    
    return nothing
end



function quiescence_differentiation!(model::ABM, cell::ActivableStemCell, event_time::Float64)
    cell.quiescence_time = draw_quiescence_time(model, 1)[1] + event_time
    cell.differentiation_time = draw_differentiation_time(model, 1)[1] + event_time
    cell.division_time = draw_division_time(model, 1)[1] + event_time

    if (cell.quiescence_time < cell.differentiation_time) & (cell.quiescence_time < cell.division_time)
        cell.next_event_time = cell.quiescence_time
        cell.next_event_type = "Quiescence"

    elseif (cell.differentiation_time < cell.quiescence_time) & (cell.differentiation_time < cell.division_time)
        cell.next_event_time = cell.differentiation_time
        cell.next_event_type = "Differentiation"
    
    elseif (cell.division_time < cell.quiescence_time) & (cell.division_time < cell.differentiation_time)
        cell.next_event_time = cell.division_time
        cell.next_event_type = "Division"

    end
    
    return nothing
end


function differentiation!(model::ABM, cell::ActivableStemCell)
    model.cell_differentiated[model.s] += 1

    cell.type = 2
    
    cell.quiescence_time = 0.0
    cell.differentiation_time = 0.0
    cell.activation_time = 0.0
    cell.division_time = 0.0

    cell.next_event_time = 0.0
    cell.next_event_type = "None"
    return nothing
end


function quiescence!(model::ABM, cell::ActivableStemCell)
 
    cell.type = 0

    cell.quiescence_time = 0.0
    cell.differentiation_time = 0.0
    cell.division_time = 0.0

    cell.activation_time = draw_activation_time(model, 1)[1] + cell.next_event_time
    cell.next_event_time = cell.activation_time
    cell.next_event_type = "Activation"

    if cell.next_event_time < model.s + 1
        agent_step!(cell, model)
    end
    
    return nothing
end


function division!(model::ABM, cell::ActivableStemCell)
    next_division_time = draw_division_time(model, 2)
    model.cell_created[model.s] += 1
    model.stem_cell_number[2][model.s] += 1
    daughters_id = [0, 0]

    for i in 1:2  
        if (cell.quiescence_time < cell.differentiation_time) & (cell.quiescence_time < next_division_time[i] + cell.next_event_time)
            daughters_id[i] = nextid(model)
            daughter = ActivableStemCell(daughters_id[i], 1, next_division_time[i] + cell.next_event_time, cell.activation_time, cell.quiescence_time, cell.differentiation_time, cell.label, cell.quiescence_time, "Quiescence")
            add_agent!(daughter, model)

        elseif (cell.differentiation_time < cell.quiescence_time) & (cell.differentiation_time < next_division_time[i] + cell.next_event_time)
            daughters_id[i] = nextid(model)
            daughter = ActivableStemCell(daughters_id[i], 1, next_division_time[i] + cell.next_event_time, cell.activation_time, cell.quiescence_time, cell.differentiation_time, cell.label, cell.differentiation_time, "Differentiation")
            add_agent!(daughter, model)
        
        elseif (next_division_time[i] + cell.next_event_time < cell.quiescence_time) & (next_division_time[i] + cell.next_event_time < cell.differentiation_time)
            daughters_id[i] = nextid(model)
            daughter = ActivableStemCell(daughters_id[i], 1, next_division_time[i] + cell.next_event_time, cell.activation_time, cell.quiescence_time, cell.differentiation_time, cell.label, next_division_time[i] + cell.next_event_time, "Division")
            add_agent!(daughter, model)
        end

        if model[daughters_id[i]].next_event_time < model.s + 1
            agent_step!(model[daughters_id[i]], model)
        end
    end
    remove_agent!(cell, model)
    
    return nothing
end

############################################################################################################
# Draw from Distributions
############################################################################################################

function draw_quiescence_time(model::ABM, n::Int)
    numbers = rand(Normal(model.quiescence_parameters[1][1],model.quiescence_parameters[1][2]), n)
    for i in 1:n
        if numbers[i] < 0
            numbers[i] = -numbers[i]
        end
    end
    return numbers
end

function draw_activation_time(model::ABM, n::Int)
    numbers = rand(Normal(model.activation_parameters[1][1], model.activation_parameters[1][2]), n)
    for i in 1:n
        if numbers[i] < 0
            numbers[i] = -numbers[i]
        end
    end
    return numbers
end

function draw_differentiation_time(model::ABM, n::Int)
    numbers = rand(Normal(model.activation_parameters[2][1], model.activation_parameters[2][2]), n)
    for i in 1:n
        if numbers[i] < 0
            numbers[i] = -numbers[i]
        end
    end
    return numbers
end

function draw_division_time(model::ABM, n::Int)
    numbers =  rand(Normal(model.division_parameters[1][1], model.division_parameters[1][2]), n)
    for i in 1:n
        if numbers[i] < 0
            numbers[i] = -numbers[i]
        end
    end
    return numbers
end


############################################################################################################
# Model running 
############################################################################################################

function agent_step!(cell::ActivableStemCell, model::ABM)
    if cell.next_event_type == "Quiescence"
        quiescence!(model, cell)
        #print("Quiescence")
    
    elseif cell.next_event_type == "Activation"
        activation!(model, cell)
    
    elseif cell.next_event_type == "Differentiation"
        differentiation!(model, cell)
    
    elseif cell.next_event_type == "Division"
        division!(model, cell)
    end
    
    return nothing
end


function model_step!(model::ABM)
    if model.s == 42
        labelling(model, 0.01, 0)
    end
    ids = collect(allids(model))
    model.stem_cell_number[1][model.s] = length(filter(id -> isequal(model[id].type, 0), ids))
    model.stem_cell_number[2][model.s] = length(filter(id -> isequal(model[id].type, 1), ids))
    model.stem_cell_number[3][model.s] = length(filter(id -> isequal(model[id].type, 2), ids))
    model.s += 1
end 

function custom_collect_stem_agent_data!(model::ABM, properties::Vector)
    dd = DataFrame()
    alla = collect(allagents(model))
    filter!(ag -> ag.type == 2, alla)
    dd[!, :step] = map(a -> model.s, alla)
    for fn in properties[2:length(properties)]
        # Extract data for each property using getfield
        dd[!, fn] = getfield.(alla, fn)
    end
    return dd
end


function custom_stem_run!(model::ABM, agent_step!::Function, model_step!::Function, n_steps::Int, properties::Vector)

    df = DataFrame()
    for i in 1:n_steps
        step!(model, agent_step!, model_step!, 1)
        dd = custom_collect_stem_agent_data!(model, properties)
        df = vcat(df, dd)
    end
    return df
end
   

function ms(model::ABM)

    ids = collect(allids(model))
        # filter all ids whose agents have `w` less than some amount
    filter!(id -> floor(model[id].next_event_time) == model.s, ids)
    return ids
end

function initialize_layersnichemodel!(n_stem_cells::Int, division_parameters::Vector, activation_parameters::Vector, quiescence_parameters::Vector, nbr_steps::Int, ms)
    space = nothing
    modelparameters_0 = LayersNicheParameters(s = 1, stem_cell_number = [[0 for _ in 1:nbr_steps], [0 for _ in 1:nbr_steps], [0 for _ in 1:nbr_steps]], activation_parameters = activation_parameters, division_parameters = division_parameters, quiescence_parameters = quiescence_parameters, cell_created = [0 for _ in 1:nbr_steps], cell_differentiated = [0 for _ in 1:nbr_steps])
    niche = StandardABM(ActivableStemCell, space; properties = modelparameters_0, scheduler = ms)
    activation_times = draw_activation_time(niche, n_stem_cells) .* rand(n_stem_cells)
    
    for i in 1:n_stem_cells
        add_agent!(ActivableStemCell, niche, 0, 0.0, activation_times[i], 0.0, 0.0, false, activation_times[i], "Activation")
    end
    
    return niche
end


############################################################################################################
# Analyses
############################################################################################################

function estimate_probability(quiescence_distribution, division_distribution, differentiation_distribution, n_samples)
    probability = Dict()
    quiescence_first = 0
    differentiation_first = 0
    division_1 = 0
    division_2 = 0
    division_3 = 0


    for _ in 1:n_samples
        # Draw samples from both distributions
        quiescence_time = rand(quiescence_distribution)
        differentiation_time = rand(differentiation_distribution)
        division_times = [rand(division_distribution)]
        push!(division_times, division_times[1] + rand(division_distribution))
        push!(division_times, division_times[2] + rand(division_distribution))

        # Check if the cells returns to quiescence
        if (quiescence_time < division_times[1]) & (quiescence_time < differentiation_time)
            quiescence_first += 1

        elseif (differentiation_time < division_times[1]) &  (differentiation_time < quiescence_time)
            differentiation_first += 1
        elseif (division_times[1] < differentiation_time) &  (division_times[1] < quiescence_time)
            division_1 += 1
            if (division_times[2] < differentiation_time) &  (division_times[2] < quiescence_time)
                division_2 += 1
                
                if (division_times[3] < differentiation_time) &  (division_times[3] < quiescence_time)
                    division_3 += 1
                end
            end
        end
        
    end

    # Estimate probability
    probability["quiescence"] = quiescence_first / n_samples
    probability["differentiation"] = differentiation_first / n_samples
    probability["1division"] = division_1 / n_samples
    probability["2division"] = division_2 / n_samples
    probability["3division"] = division_3 / n_samples
    return probability
end


function labelling(model::ABM, frequency::Float64, type::Int)
    ids = collect(allids(model))
    filter!(id -> model[id].type == type, ids)
    #print(ids)
    labeled_cells = rand(ids, Int(floor(frequency*length(ids))))
    for i in labeled_cells
        model[i].label = true 
    end
    return labeled_cells
end

#data collection
type0(a) = a.type==0
type1(a) = a.type==1
type2(a) = a.type==2
labeled2(a) = (a.label==true && a.type==2)
labeled1(a) = (a.label==true && a.type==1)
labeled0(a) = (a.label==true && a.type==0)


#Adam mac Lean M3 model 2014 Paper

#rho = 1/110 #growth rate 
#K = 17000   #Carrying capacity
#delta = 0.1 # "death rate" actually its more a rate of "leaving the stem cell state wether through death or differentiation 
#phi = x, rho, K -> (rho + 1)*K/(K+(rho*x))
#function phi(rho, K, x)
 #   return (rho + 1)*K/(K+(rho*x))
#end
#p = [rho, K, delta]
#function log_growth!(du,u,p,t) 
#    phi = (p[1] + 1)*p[2]/(p[2]+(p[1]*u[1]))
#    du[1] = (phi - 1)*p[3]*u[1]
#    du[2] = p[3]*u[1]
#end

#prob = ODEProblem(log_growth!, [1580.0; 0.0], (0.0, 15000.0), p = (0.1, 1700, 0.01))
#sol = solve(prob)