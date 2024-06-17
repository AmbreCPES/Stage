parameters_file = ""
output_file = ""
n_steps = 10
count_file = ""
histo_file = ""
plt_file = ""
test_mean=50

print("Arguments provided to execute_run_ABM.jl [n_steps, parameters_file, output_file, count_file, histo_file, plt_file]:", ARGS)
if length(ARGS) == 7
    n_steps=parse(Int64, ARGS[1])
    parameters_file = ARGS[2]
    output_file = ARGS[3]
    count_file = ARGS[4]
    histo_file = ARGS[5]
    plt_file = ARGS[6]
    test_mean = parse(Int64, ARGS[7])
end


include("stem_cell.jl")
using CSV, DelimitedFiles

if parameters_file != ""

include(parameters_file)
    
if @isdefined(division) == false
    error("\n transition_matrix_0 is not defined in parameters_file")
end

if @isdefined(activation) == false
    error("\n ttd is not defined in parameters_file. Impossible to initialize TypeParameters")
end

if @isdefined(differentiation) == false
    error("\n ttnd is not defined in parameters_file. Impossible to initialize TypeParameters")
end

if @isdefined(quiescence) == false
    error("\n collection_t0 is not defined in parameters_file. No cells defined at t = 0")
end


#Data collection
#type0 quiescent HSCs
#type1 active HSCs
#type2 differentiated cells

layersniche = initialize_layersnichemodel!(17000, [(division[1], division[2])], [(test_mean, activation[2]), (differentiation[1], differentiation[2])], [(quiescence[2], quiescence[2])], n_steps, ms)
data = custom_stem_run!(layersniche, agent_step!, model_step!, n_steps, [:step, :id, :type, :label]) # the first symbol has to be :step 
println(layersniche.activation_parameters)

CSV.write(output_file, data)
matrix_count = hcat(layersniche.stem_cell_number[1], layersniche.stem_cell_number[2], layersniche.stem_cell_number[3])
#count_type0, count_type1, count_type2 = layersniche.stem_cell_number[1], layersniche.stem_cell_number[2], layersniche.stem_cell_number[3]
io = open(count_file, "w") do io
    writedlm(io, matrix_count)
end


if !isfile(histo_file)
    numbers = draw_quiescence_time(layersniche, 1000)
    numbers1 = draw_activation_time(layersniche, 1000)
    numbers2 =  draw_division_time(layersniche, 1000)
    numbers3 = draw_differentiation_time(layersniche, 1000)
    histogram(numbers, bins=20, xlabel="Value", ylabel="Frequency", label="quiescence time", title="Histogram of Time distribution")
    histogram!(numbers1, bins=20, label="acivation time")
    histogram!(numbers2, bins=20, label="division time")
    histogram!(numbers3, bins=20, label="differetiation time")
    savefig(histo_file)
end

#plot(1:n_steps, layersniche.stem_cell_number[1] + layersniche.stem_cell_number[2], label = "Total HSCs", title = "Evolution of HSCs in the niche", xlabel="Time(in days)", ylabel="Number of cells", show = true)
#plot!(1:n_steps, layersniche.stem_cell_number[2] , label = "Activated HSC")
#plot!(1:n_steps, layersniche.stem_cell_number[3] , label = "Differentiated cells")
#ylims!(0,25000)
#savefig(plt_file)

end