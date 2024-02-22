using Plots # you should really use Makie but that's for later
using Random
using StatsBase

# DISCLAIMER: very unoptimised at the moment, happy to chat about optimising later

# simulate one time step of the branching process
# offspring_func : function that randomly generates number of offspring for one individual
function step_bp(pop_size, offspring_func, death_rate_t)
	new_pop_size = 0
	for _ in 1:pop_size
        if death_func(death_rate_t) == 1
		    new_pop_size += offspring_func(pop_size)
        end
	end
	return new_pop_size
end


#Okay Now I don't want the population to stop growing because of a carrying capacity but because cell death rate increase exponentially with time
function death_func(death_rate)
	death_prob = [death_rate, 1 - death_rate] # probability of dying, and surviving
	return sample(0:length(death_prob)-1, ProbabilityWeights(death_prob))
end

function death_rate(t, a=0.009, c=0)
	# death rate follows a linear pattern dependant of time : at+c
	return a*t+c
end

function new_offspring_func(pop_size)
	offspring_prob = [0.9, 0.1] # probability of staying put and dividing
	return sample(1:length(offspring_prob), ProbabilityWeights(offspring_prob))
end

function simulate_bp(steps, offspring_func)
	bp = zeros(steps)
	bp[1] = 1
	for t in 2:steps
		death_rate_t = death_rate(t)
		bp[t] = step_bp(bp[t-1], offspring_func, death_rate_t)#ici on regarde la population à t-1 donc ici on fix pop_size=1 à t=0
        if bp[t] == 0
            print("Everyone died :( in ", t, " time steps")
            break
        end
	end
	return bp
end



# simulate 10 independent branching processes
bps = [simulate_bp(100, const_offspring_func) for _ in 1:10];
plot(bps, label=nothing) # these explode or go extinct


#Actually this seems trivial but make sens for what we are trying to model
#When we look at HSC transplantation often you have an initial increase but then graft fails
#because the number of "young MPP" is not growing fast enough and so eventually differientated cell population crash after a few weeks

#So now what if we tried to implement self renewal in this process
# death rate is dependant of cells' age and not of the system's age