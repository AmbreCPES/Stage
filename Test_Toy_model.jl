include("Function_Toy_model.jl")





cell0 = [[1], [2]]
life = initialize_model(cell0, ModelParameters())

activation_rate = life.activation_rate
activation = rand(Binomial(life.pool, activation_rate), 1)

for _ in 1:activation[1] 
    hc = HematopoeiticCell(nextid(life), 1, 0, 0, 0, [], round(shifted_exponential_sampling(life.division_param[1], life.division_param[2])))
    add_agent!(hc, life)
end

cells_creation!(life, 100)

life