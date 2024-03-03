#Univariate Branching process using Agents.jl
using Agents 
using StatsBase
using Plots
using BenchmarkTools
using Statistics
using DataFrames

@agent HSC NoSpaceAgent begin
    quiescence::Bool # whether the agent is likely to divide True = quiescence
    age::Int
end

@agent hemato_cell NoSpaceAgent begin
    age::Int #Age of the cell since she left HSC in number of time step
    nbr_division::Int #nbr of divisions since the Cell left HSC
end


function initialize_model(nbr_HSC,nbr_hemato)
    properties = Dict(:Time => 0)
    space = nothing 
    life = ABM(Union{HSC, hemato_cell}, space, scheduler = Schedulers.ByType(false, true, Union{HSC, hemato_cell}); properties)

    first_HSC = HSC.(1:nbr_HSC, trues(nbr_HSC), zeros(nbr_HSC))
    first_hemato = hemato_cell.(1:nbr_hemato, zeros(nbr_hemato), zeros(nbr_hemato))
    for cell in first_HSC 
        add_agent!(cell, life)
    end
    for cell in first_hemato 
        add_agent!(cell, life)
    end
    return life
end


function division!(cell::hemato_cell, model; division_rate = 0.3)
    
    if sample(0:1, ProbabilityWeights([(1 - division_rate), division_rate])) == 1
        divisions = cell.nbr_division + 1

        add_agent!(hemato_cell, model, cell.age, divisions)
        add_agent!(hemato_cell, model, cell.age, divisions)
        remove_agent!(cell, model)
    end
end

function division!(cell::HSC, model; division_rate = 0.5)
    if sample(0:1, ProbabilityWeights([(1 - division_rate), division_rate])) == 1
        
        add_agent!(HSC, model, cell.quiescence, cell.age)
        add_agent!(HSC, model, cell.quiescence, cell.age)
        remove_agent!(cell, model)
    end
end

#Les fonctions divisions fonctionnent pour HSC et hemato_cell

function cell_activation!(cell::HSC, activation_rate)
    if sample(0:1, ProbabilityWeights([1 - activation_rate, activation_rate])) == 1
        cell.quiescence = false  
    end
end
#fonctionne

function cell_transition!(cell::HSC, model, transition_rate)
    if sample(0:1, ProbabilityWeights([1 - transition_rate, transition_rate])) == 1
        add_agent!(hemato_cell, model, 0, 0)
        remove_agent!(cell, model)
    end
end

#fonctionne

function cells_creation(model, creation_rate = 0.5)
    create = sample(0:1, ProbabilityWeights([1 - creation_rate, creation_rate]))
    if create == 1
        add_agent!(hemato_cell, model, 0, 0)
    end
end

function death_func(age, a=0.009, c=0)
    death_rate = a * age + c
	death_prob = [death_rate, 1 - death_rate]
	return sample(0:1, ProbabilityWeights(death_prob))
end

function life_step!(cell::HSC, model)
    
    if cell.quiescence == true
        cell.age += 1 
        cell_activation!(cell, 0.02)
    else
        cell_transition!(cell, model, 1)
    end
    if typeof(cell) == HSC
        division!(cell, model; division_rate = 0.5)
    end
end
# Fonctionne

function life_step!(cell::hemato_cell, model)
    
    if death_func(cell.age) == 1 #in case of survival
        division_rate = 0.3
        cell.age += 1
        division!(cell, model; division_rate)
        #print("survival")
    else
        remove_agent!(cell, model)     
    end
    
end

function model_step!(model)
    model.Time += 1
    cells_creation(model)
end

hsc(a) = typeof(a) == HSC
hemato(a) = typeof(a) == hemato_cell
adata = [(hsc, count), (hemato, count)]
mdata = [:Time]

n_steps= 500
df = DataFrame(step =[0], count_hsc =[0], count_hemato = [0])
for _ in 1:10
    life = initialize_model(0,1)
    data,_ = run!(life, life_step!, model_step!, n_steps; adata, mdata)
    append!(df, data)
 end
 pretty = [df.count_hemato[(i-1)*n_steps+i+1:((i-1)*n_steps+n_steps+1+i), : ] for i in 1:10]
 plot(1:(n_steps+1), pretty, labels =["hemato"])

#plot(data.step, [data.count_hemato, data.count_hsc, data.count_hemato + data.count_hsc], labels =["hemato" "HSC" "total"])