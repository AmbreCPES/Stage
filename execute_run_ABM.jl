parameters_file = ""
output_file = ""
output_death_file = ""

Base.@kwdef struct TypeParameters
    # each of the properties has for length the number of type
    ttd::Vector{Tuple} = [(10000.0, 0.2), (330.0, 0.2), (70.0, 0.2), (60.0, 0.2), (1.1, 0.2)] #Time to death
    ttnd::Vector{Tuple} = [(110, 0.1), (1/0.042, 0.2), (0.25, 0.1), (86.66, 0.18), (10.9, 0.23)] # Time to next division
end

Base.@kwdef struct TransitionDistribution
    
    ttnt::Vector{Vector{Float64}} = []
end


include("Initialization.jl")

if parameters_file == ""
    n_steps = 100
    typeparameters_0 = TypeParameters()

    #From data of supplemental material of Fundamental properties of unperturbed haematopoiesis from stem cells in vivo
    #if we know n_HSC we can define relative compartment size. 
    n_HSC = 1.0
    n_STHSC = 2.9
    n_MPP = 9.0
    n_CMP = 39.0
    n_CLP = 13.0

    #Differentiation rate:
    a_HSC_STHSC = 0.009
    a_STHSC_MPP = 0.045
    a_MPP_CMP = 3.992
    a_MPP_CLP = 0.022

    # Transition Matrix = [HSC, ST-HSC, MPP, CLP, CMP, ; ]

    transition_matrix_0 = [(1 - a_HSC_STHSC/n_HSC) (a_HSC_STHSC/n_HSC) 0.0 0.0 0.0 ;   #HSC
                        0.0 (1 - a_STHSC_MPP/n_STHSC) (a_STHSC_MPP/n_STHSC) 0.0 0.0 ;       #ST-HSC
                        0.0 0.0 (1 - a_MPP_CLP/n_MPP - a_MPP_CMP/n_MPP) a_MPP_CLP/n_MPP a_MPP_CMP/n_MPP ;  #MPP
                        0.0 0.0 0.0 1.0 0.0 ;  #CLP
                        0.0 0.0 0.0 0.0 1.0 ;          #CMP
                        ]   
                        
    n_tot = 30
    collection_t0 = [(1,10), (2,20)]

else
    include(parameters_file)
    if @isdefined(transition_matrix_0) == false
        error("\n transition_matrix_0 is not defined in parameters_file")
    end

    if @isdefined(ttd_0) == false
        error("\n ttd is not defined in parameters_file. Impossible to initialize TypeParameters")
    end

    if @isdefined(ttnd_0) == false
        error("\n ttnd is not defined in parameters_file. Impossible to initialize TypeParameters")
    end

    if @isdefined(collection_t0) == false
        error("\n collection_t0 is not defined in parameters_file. No cells defined at t = 0")
    end

    if @isdefined(n_steps) == false
        n_steps = 365
        print("\n n_steps undefined in parameters_file. Value n_steps=365 was attributed by default")
    end
end

n_tot = sum([collection_t0[x][2] for x in eachindex(collection_t0)])
typeparameters_0 = TypeParameters(ttd = ttd_0, ttnd = ttnd_0)

if size(transition_matrix_0)[1] == length(typeparameters_0.ttd) == length(typeparameters_0.ttnd)

    modelparameters_0 = initialize_modelparameters(n_steps, size(transition_matrix_0)[1], output_death_file, transition_matrix_0, typeparameters_0)
    life = initialize_model(collection_t0, transition_matrix_0, modelparameters_0, typeparameters_0, n_tot, ms)

    data = custom_run!(life, life_step!, model_step!, n)
    CSV.write(output_file, data)
else
    error("Problem with the size of transition_matrix_0, ttd_0 or ttnd_0 defined in parameters_file")
end

