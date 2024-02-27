#Univariate Branching process using Agents.jl
using Agents 
using StatsBase

@agent HSC NoSpaceAgent begin
    quiescence::Bool # whether the agent is likely to divide True = quiescence
    g::Int
end

@agent hemato_cell NoSpaceAgent begin
    age::Int #Age of the cell since she left HSC in number of time step
    nbr_division::Int #nbr of divisions since the Cell left HSC
    g::Int
end


function initialize_model(nbr_HSC)
    properties = Dict(:Time => 0)
    space = nothing 
    life = ABM(Union{HSC, hemato_cell}, space, scheduler = Schedulers.ByType(false, true, Union{HSC, hemato_cell}); properties)

    first_agents = HSC.(1:nbr_HSC, trues(nbr_HSC), ones(nbr_HSC))
    for cell in first_agents
        add_agent!(cell, life)
    end
    #print(nagents(life))
    return life
end


function division!(cell::hemato_cell, model; division_rate = 0.3)
    
    if sample(0:1, ProbabilityWeights([(1 - division_rate), division_rate])) == 1
        divisions = cell.nbr_division + 1

        add_agent!(hemato_cell, model, cell.age, divisions, 1)
        add_agent!(hemato_cell, model, cell.age, divisions, 1)
        remove_agent!(cell, model)
    end
end

function division!(cell::HSC, model; division_rate = 0.3)
    if sample(0:1, ProbabilityWeights([(1 - division_rate), division_rate])) == 1
        
        add_agent!(HSC, model, cell.quiescence, 0)
        add_agent!(HSC, model, cell.quiescence, 0)
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
        add_agent!(hemato_cell, model, 0, 0, 1)
        remove_agent!(cell, model)
    end
end

#fonctionne


function death_func(age, a=0.009, c=0)
    death_rate = a * age + c
	death_prob = [death_rate, 1 - death_rate]
    #print(death_prob) # probability of dying, and surviving
	return sample(0:1, ProbabilityWeights(death_prob))
end

function life_step!(cell::HSC, model)
    
    if cell.quiescence == true 
        cell_activation!(cell, 0.5)
        division!(cell, model; division_rate = 0.3)
    else 
        cell_transition!(cell, model, 1)
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


n_steps = 10
life = initialize_model(1)
adata = [:g]
#hsc(a) = a isa HSC
#wolf(a) = a isa hemato_cell
#adata = [(HSC, count), (wolf, count)]
data, _ = run!(life, life_step!, n_steps; adata)

