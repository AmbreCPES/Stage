#Multivariate Branching process using Agents.jl
using Agents 
using StatsBase
using Plots
using BenchmarkTools
using Statistics
using DataFrames
using DelimitedFiles

# Contradiction between Agents.jl latest doc and julia peroformace tip doc about the most efficient way to define agents

@agent hemato_cell NoSpaceAgent begin
    type::Int # A number between 1 and 20 
    age::Int #Age of the cell since she left HSC in number of time step
    nbr_division::Int #nbr of divisions since the Cell left HSC
    lineage::Array{Int64, 1}
end

#More efficient to only have one type of agents BUT ! a trade-off with memory allocations if I have more parameters 
#Must check run time and memory allocation to see what's better


#Model Initialization:

function initialize_model(t0_cells, Parameters) # We expect t0_cells to be an a vector of arrays , each array is a cell with given parameters (type, age, nbr_divisions)
    space = nothing 
    life = ABM(hemato_cell, space; properties = Parameters)
    first_cells = hemato_cell.(1:length(t0_cells),map(x -> x[1], cell0), map(x -> x[2], cell0),map(x -> x[3], cell0),)
    
    for cell in  first_cells
        add_agent!(cell, life)
    end
    return life
end

function cells_creation(model)
    creation_rate = model.creation_rate
    create = sample(0:1, ProbabilityWeights([1 - creation_rate, creation_rate]))
    if create == 1
        add_agent!(hemato_cell, model, 0, 0)
    end
end

function death_func(age, a = 0.009, c = 0)
    death_rate = a * age + c
	death_prob = [death_rate, 1 - death_rate]
	return sample(0:1, ProbabilityWeights(death_prob))
end

function division_func(age, a = -0.009, c = 0)
    division_rate = a * age + c
	division_prob = [division_rate, 1 - division_rate]
	return sample(0:1, ProbabilityWeights(division_prob))
end

function division!(cell::hemato_cell, model, division_rate)
    if sample(0:1, ProbabilityWeights([(1 - division_rate), division_rate])) == 1
        divisions = cell.nbr_division + 1

        add_agent!(hemato_cell, model, cell.age, divisions)
        add_agent!(hemato_cell, model, cell.age, divisions)
        remove_agent!(cell, model)
    end
end

function life_step!(cell::hemato_cell, model)
    
    if death_func(cell.age) == 1 #in case of survival
        cell.age += 1
        division!(cell, model, division_rate)
    else
        remove_agent!(cell, model)     
    end
    
end


#Testing:
#Parameters definition and initialization:
Base.@kwdef mutable struct Parameters
    creation_rate::Float64 = 1.0
    deaths::DataFrame
end
cell0 = [[1,0,0],[2,0,0]]
file = 
#initialize_model(cell0, Parameters)