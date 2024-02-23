using Plots # you should really use Makie but that's for later
using Random
using StatsBase

function population_init(nbr_cells)
    cells_id = collect(1:nbr_cells) #array of integers from 1 to nbr_cells
    cells_age = zeros(nbr_cells)
    return [cells_id, cells_age]
end

function new_offspring_func()
	offspring_prob = [0.9, 0.1] # probability of staying put and dividing
	return sample(1:length(offspring_prob), ProbabilityWeights(offspring_prob))
end


function cells_creation(creation_rate)
    sample(0:1, ProbabilityWeights([1 - creation_rate, creation_rate]))
end

function death_func(age, a=0.009, c=0)
    death_rate = a * age + c
	death_prob = [death_rate, 1 - death_rate] # probability of dying, and surviving
	return sample(0:length(death_prob)-1, ProbabilityWeights(death_prob))
end

function step_bp_2(tot_cells, cells_id, cells_age, offspring_func, creation_rate)
    new_cells = []
    dead_cells = []
    pop_size = length(cells_id)
	for i in 1:pop_size
        age = cells_age[i]
        if death_func(age) == 1
		    if offspring_func(pop_size) == 1
                cells_age[i] += 1
            else
                append!(new_cells,[age + 1, age + 1])
                append!(dead_cells, i)
            end  
        else
            append!(dead_cells, i)
        end
	end

    if cells_creation(creation_rate) == 1
        append!(new_cells, 0)
    end

    #print("dead_cells = ", dead_cells, " new_cells = ", new_cells)
    deleteat!(cells_id, dead_cells) 
    deleteat!(cells_age, dead_cells)

    append!(cells_id, (tot_cells + 1):(tot_cells + length(new_cells)))
    append!(cells_age, new_cells)
    #print("cells_id = ", cells_id, " cells_age = ", cells_age)
	return cells_id, cells_age, length(new_cells), length(dead_cells)
end



function simulate_bp(steps, offspring_func, nbr_cells_0, creation_rate)
	
    bp = zeros(steps)
    nbr_dead = 0
    tot_cells = nbr_cells_0
	bp[1] = nbr_cells_0

    cells_id, cells_age = population_init(nbr_cells_0)

	for t in 2:steps
		cells_id, cells_age, new_cells, dead_cells = step_bp_2(tot_cells, cells_id, cells_age, offspring_func, creation_rate)#ici on regarde la population à t-1 donc ici on fix pop_size=1 à t=0
        bp[t] = length(cells_id)
        tot_cells += new_cells
        nbr_dead += dead_cells
        if bp[t] == 0
            print("Everyone died :( in ", t, " time steps")
            break
        end
	end
	return bp
end


#step_bp_2(1, [1], [0], new_offspring_func, 0.5)

bps = [simulate_bp(500, new_offspring_func, 1, 0.5) for _ in 1:10];
plot(bps, label=nothing) # these explode or go extinct