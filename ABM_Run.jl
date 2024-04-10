
include("Initialization.jl")

#= Parameter setting from the article hematopoeisis in number, we are going to use the mice data
Types and their characteristics:

Total number of hematopoeitic cell in a mouse = 1.7 *10^10 ( 1.6*10^10 in blood, 4.5*10^8 in bone marrow, 2.1*10^8 in thymus, 7.6*10^7 in lymph node, 2.1*10^8 in the spleen)
Total cell number across different compartment in a healthy mice: 10^4 HSC, 10^5 MPP, 10^6 RPP, 10^8 Granulomyeloids, 10^9 Platelets, 10^10 Erythrocytes, 10^8 Lymphoids

MPP: type 1: 

Information found in supplementary material 1:

Lifespan = 1/Death rate = t1/2/ln(2)
half life = time for half to die 
B Cell: Half-life = 24,3 days, lifespan = 35,00 days, death_rate = 0.03
Granulomyeloids: half-life = 0.84, lifespan = 1.21 days, death_rate = 0.83
Platelets: Half-life = 2.5 days, lifespan = 3.60 days, death_rate = 0.28
Red blood cells: half-life = 24.3 days, lifespan = 35 days, death_rate = 0.03

Parameters settings found in the article "Fundamental properties of unperturbed haematopoiesis from stem cells in vivo"
On average, per day, 1 out of 110 HSCs differentiates into an ST-HSC, and 1 out of 22 ST-HSCs differentiates into an 
MPP. At the MPP stage, considered a lymphoid–myeloid bifurcation point, we estimate that per day 1 out of 46 MPPs generates a CLP, 
while 1 MPP generates 4 CMPs. Given that the cell numbers also increase from HSC to ST-HSC and MPP, the
efflux of cells exceeded influx in all of these compartments. To maintain compartment size, this flux difference is balanced by net proliferation
=#
Base.@kwdef struct TypeParameters1
    # each of the properties has for length the number of type
    ttd::Vector{Tuple} = [(μ_for_mean(86.66, 0.18), 0.18), (μ_for_mean(86.66, 0.18), 0.18)] #Time to death
    ttnd::Vector{Tuple} = [(μ_for_mean(10.9, 0.23), 0.23), (μ_for_mean(86.66, 0.18), 0.18)] # Time to next division
    dd::Vector{Tuple} = [(μ_for_mean(53.79, 0.21), 0.21), (μ_for_mean(86.66, 0.18), 0.18)] #division destiny
end

# Transition Matrix = [MPP4 CLP Dendritic_cell NK B_Cell T_Cell; ]

transition_matrix_0 = [0.9 0.06 0.01 0.01 0.01 0.01;
                      0.0 0.9 0.025 0.025 0.025 0.025;
                      0.0 0.0 1.0 0.0 0.0 0.0;
                      0.0 0.0 0.0 1.0 0.0 0.0;
                      0.0 0.0 0.0 0.0 1.0 0.0;
                      0.0 0.0 0.0 0.0 0.0 1.0]


#so here n states = 2                  
typeparameters_0 = TypeParameters()
n_states = length(typeparameters_0.ttd)

n_tot = 30
collection_t0 = [(1, 14), (2, 16)]

type_distribution_0 = initialize_distributions(typeparameters_0, n_states)

cell_collection_0 = initialize_collection(type_distribution_0, collection_t0, n_tot)

death_file = "C:\\Users\\ambre\\Documents\\ENS\\stage_M1\\code_1\\test_file"
modelparameters_0 = initialize_modelparameters(6, death_file, transition_matrix_0, typeparameters_0)

life = initialize_model(collection_t0, modelparameters_0, typeparameters_0, n_tot)

division!(life, life[1], 2)

