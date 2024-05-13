using DataFrames, CSV

##########################################################################################################################################
# Plotting data as graph
##########################################################################################################################################


"works only with the old way data were stored (all lineage of each cell as stored)"
function cell_count(data::DataFrame, death_file::String, cell_type::Vector{Int}, last_time_step::Int, nbr_start_cell::Int)
    deaths = DataFrame(CSV.File(death_file))
    cells_data = deepcopy(data)
    cells_data = vcat(cells_data, deepcopy(deaths))
    sort!(cells_data, :id, rev = true)
    last_id = cells_data[1,:id]
    cells_number = Dict{Int, Vector{Int}}(type => zeros(last_time_step) for type in cell_type)
    cell_to_count = [string(i)*"_"*string(type) for i in 1:last_id for type in cell_type]
    
    for cell in 1:nrow(cells_data)
        
        if cells_data[cell, :id] in data[!, :id]
            cells_data[cell, :lineage] = vcat(cells_data[cell, :lineage], [cells_data[cell, :id],last_time_step, cells_data[cell, :type]])
        else
            cells_data[cell, :lineage] = parse.(Int, split(chop(cells_data[cell, :lineage], head = 1, tail = 1), ","))
        end
    
        lineage = cells_data[cell,:lineage]
    
        type = lineage[length(lineage)]
        time = lineage[length(lineage) - 1]
        id = lineage[length(lineage) - 2]
    
            if length(lineage) >= 6
                for i in range(start = (length(lineage) - 5), stop = 0, step = -3)
                    
                    if string(id)*"_"*string(type) in cell_to_count
                        cells_number[type][(lineage[i + 1] + 1):time] = cells_number[type][(lineage[i + 1] + 1):time] .+ 1
                        filter!(!=(string(id)*"_"*string(type)), cell_to_count)
                    else
                        break
                    end
                    
                    type = lineage[i + 2]
                    time = lineage[i + 1]
                    id = lineage[i]
                end
                if (string(id)*"_"*string(type) in cell_to_count) & (id <= nbr_start_cell)
                    cells_number[lineage[3]][1:time] = cells_number[lineage[3]][1:time] .+ 1
                    filter!(!=(string(id)*"_"*string(type)), cell_to_count)
                end
            
            elseif length(lineage) == 3
                if string(id)*"_"*string(type) in cell_to_count
                    
                    cells_number[lineage[3]][1:time] = cells_number[lineage[3]][1:time] .+ 1
                    filter!(!=(string(id)*"_"*string(type)), cell_to_count)
                end
            end
        
    end
    
    
    return cells_number #a dict with key the type and associated a vector with the number of cell of this type at each time step
end

"""
cell_count2(data::DataFrame, death_file::String, cell_type::Vector{Int}, last_time_step::Int, nbr_start_cell::Int)

Special case of precision_cell_count where precision = 1

"""

function cell_count2(data::DataFrame, death_file::String, cell_type::Vector{Int}, last_time_step::Int, nbr_start_cell::Int)
    deaths = DataFrame(CSV.File(death_file))
    cells_data = deepcopy(data)
    cells_data = vcat(cells_data, deepcopy(deaths))
    sort!(cells_data, :id)
    cells_number = Dict{Int, Vector{Int}}(type => zeros(last_time_step + 1) for type in cell_type)
    cell_to_count = [string(i)*"_"*string(type) for i in 1:nbr_start_cell for type in cell_type]
    
    for cell in 1:nrow(cells_data)
        
        if cells_data[cell, :id] in data[!, :id]
            cells_data[cell, :lineage] = vcat(cells_data[cell, :lineage], [cells_data[cell, :id], last_time_step, cells_data[cell, :type]])
        else
            cells_data[cell, :lineage] = parse.(Float64, split(chop(cells_data[cell, :lineage], head = 1, tail = 1), ","))
        end
    
        lineage = cells_data[cell,:lineage]
    
        type = Int(lineage[length(lineage)])
        time = Int(floor(lineage[length(lineage) - 1]))
        id = Int(lineage[length(lineage) - 2])
    
            if length(lineage) >= 6
                for i in range(start = (length(lineage) - 5), stop = 0, step = -3)
                        
                    cells_number[type][(Int(floor(lineage[i + 1])) + 1):(time+1)] = cells_number[type][(Int(floor(lineage[i + 1])) + 1):(time+1)] .+ 1
                    
                    type = Int(lineage[i + 2])
                    time = Int(Int(floor(lineage[i + 1])))
                    id = Int(lineage[i])
                end

                if (string(id)*"_"*string(type) in cell_to_count)
                    cells_number[type][1:time] = cells_number[type][1:time] .+ 1
                    filter!(!=(string(id)*"_"*string(type)), cell_to_count)
                end
            
            elseif length(lineage) == 3
                if string(id)*"_"*string(type) in cell_to_count
                    
                    cells_number[type][1:time] = cells_number[type][1:time] .+ 1
                    filter!(!=(string(id)*"_"*string(type)), cell_to_count)
                end
            end
        
    end
    
    
    return cells_number #a dict with key the type and associated a vector with the number of cell of this type at each time step
end

"""
function precision_cell_count(
    data::DataFrame : the output of custom_run!
    death_file::String :  file in which dead cells information where registered during custom_run!
    cell_type::Vector{Int} : a vector of Int containing all existing cell type. 
    last_time_step::Int : time steps at which the models stops, either equal to model.n_steps or to the first time step where all cells are dead 
    nbr_start_cell::Int : we indicate the number of cell we initialize in the model, we assume that if 100 cells are present at time 0 they are attributed the ids 1 to 100
    precision::Float64 : a number equal to 1, 0.1, 0.001, ... indicated the precision of the count. If precision = 1, we get the number of cell per day, if 0.1 we get the number of cell per 0.1 day
    )

This functions takes the output of custom run, a Dataframe and a death_file, and returns a dict with keys the cells type and associated to them a vector.
This vector's length is the number of time step in the run multiplied by 1/precision. For instance, if precision is 0.1 and the model ran for 365 the vectors will be 3660 long.
It's not 3650 because we count also between 0 and 1 the model starting at day 0. 
"""
function precision_cell_count(file, death_file::String, cell_type::Vector{Int}, last_time_step::Int, nbr_start_cell::Int, precision::Float64)
    
    if isa(file, String)
        data = DataFrame(CSV.File(file))
    end

    cells_data = DataFrame(id = Int64[], type = Int64[], lineage = Any[])
    cells_data = vcat(cells_data, deepcopy(data))
    deaths = DataFrame(CSV.File(death_file)) 
    cells_data = vcat(cells_data, deepcopy(deaths))
    sort!(cells_data, :id)

    cells_number = Dict{Int, Vector{Int}}(type => zeros(length(range(start = 0, stop = last_time_step + 1, step = precision))) for type in cell_type)
    cell_to_count = [i for i in 1:nbr_start_cell]
    #print(isa(cells_data[cell, :lineage],String15))


    for cell in 1:nrow(cells_data)
        
        if cells_data[cell, :id] in data[!, :id]

            if (isa(cells_data[cell, :lineage],Vector) == false) & (cells_data[cell, :lineage] != "Float64[]")
                cells_data[cell, :lineage] = parse.(Float64, split(chop(cells_data[cell, :lineage], head = 1, tail = 1), ","))
                cells_data[cell, :lineage] = vcat(cells_data[cell, :lineage], [cells_data[cell, :id], last_time_step, cells_data[cell, :type]])

            elseif (cells_data[cell, :lineage] == "Float64[]")
                cells_data[cell, :lineage] = [cells_data[cell, :id], last_time_step, cells_data[cell, :type]]
            end
            
        elseif cells_data[cell, :id] in deaths[!, :id]
            
            cells_data[cell, :lineage] = parse.(Float64, split(chop(cells_data[cell, :lineage], head = 1, tail = 1), ","))
        else
            print("\n mais d'ou vient cette cellule ?: ", cells_data[cell, :id] )
        end

        lineage = cells_data[cell,:lineage]
        
        type = Int(lineage[length(lineage)])
        time = Int(floor(lineage[length(lineage) - 1]*(1/precision))+(1/precision))
        id = Int(lineage[length(lineage) - 2])

        if length(lineage) >= 6
            for i in range(start = (length(lineage) - 5), stop = 0, step = -3)
                        
                cells_number[type][(Int(floor(lineage[i + 1]*(1/precision)) + (1/precision))):(time-1)] = cells_number[type][(Int(floor(lineage[i + 1]*(1/precision)) + 1/precision)):(time-1)] .+ 1
                type = Int(lineage[i + 2])
                time = Int(floor(lineage[i + 1]*(1/precision))+Int(1/precision))
                id = Int(lineage[i])
            end

            if id in cell_to_count
                cells_number[type][1:(time-1)] = cells_number[type][1:(time-1)] .+ 1
                filter!(!=(id), cell_to_count)
            end
                
        elseif length(lineage) == 3
            if id in cell_to_count
                cells_number[type][1:(time)] = cells_number[type][1:(time)] .+ 1
                filter!(!=(id), cell_to_count)
            end
        end    
    end

    return cells_number
end

function precision_cell_count2(file, death_file::String, cell_type::Vector{Int}, last_time_step::Int, nbr_start_cell::Int, precision::Float64)
    
    if isa(file, String)
        data = DataFrame(CSV.File(file))
    end

    cells_data = DataFrame(id = Int64[], type = Int64[], lineage = Any[])
    cells_data = vcat(cells_data, deepcopy(data))
    deaths = DataFrame(CSV.File(death_file)) 
    cells_data = vcat(cells_data, deepcopy(deaths))
    sort!(cells_data, :id)

    cells_number = Dict{Int, Vector{Int}}(type => zeros(length(range(start = 0, stop = last_time_step, step = precision))) for type in cell_type)
    cell_to_count = [i for i in 1:nbr_start_cell]
    #print(isa(cells_data[cell, :lineage],String15))


    for cell in 1:nrow(cells_data)
        
        if typeof(cells_data[cell, :lineage]) != Vector
            cells_data[cell, :lineage] = parse.(Float64, split(chop(cells_data[cell, :lineage], head = 1, tail = 1), ","))
        else
            print("\n mais d'ou vient cette cellule ?: ", cells_data[cell, :id] )
        end

        lineage = cells_data[cell,:lineage]
        
        type = Int(lineage[length(lineage)])
        time = Int(floor(lineage[length(lineage) - 1]*(1/precision))+(1/precision)) + 1
        id = Int(lineage[length(lineage) - 2])

        if length(lineage) >= 6
            for i in range(start = (length(lineage) - 5), stop = 0, step = -3)
                        
                cells_number[type][(Int(floor(lineage[i + 1]*(1/precision)) + (1/precision))):(time-1)] = cells_number[type][(Int(floor(lineage[i + 1]*(1/precision)) + 1/precision)):(time-1)] .+ 1
                type = Int(lineage[i + 2])
                time = Int(floor(lineage[i + 1]*(1/precision))+Int(1/precision))
                id = Int(lineage[i])
            end

            if id in cell_to_count
                cells_number[type][1:(time-1)] = cells_number[type][1:(time-1)] .+ 1
                filter!(!=(id), cell_to_count)
            end
                
        elseif length(lineage) == 3
            if id in cell_to_count
                cells_number[type][1:(time-1)] = cells_number[type][1:(time-1)] .+ 1
                filter!(!=(id), cell_to_count)
            end
        end    
    end

    return cells_number
end

###############################################################
#                LINEAGE NEWICK FORMATING
###############################################################

function find_root(lineage::Vector{Vector{Float64}}, root_id::Float64; cell_count::Int = 1) # une fonction qui trouve dans une lste de liste les vecteur commencant par un nombre donnée
    cell_found = []
    for x in eachindex(lineage)
        if lineage[x][1] == root_id
            push!(cell_found, x)
            if length(cell_found) == cell_count
                return cell_found
            end
        end
    end
    return cell_found
end

function find_cell(cells_data::DataFrame, cell_id)
    return findfirst(x -> x == Int(cell_id), cells_data.id)
end

function daughters_id(lineage::Vector{Vector{Float64}}, mother_row::Int, pos::Int)
    time_div = lineage[mother_row][pos + 1]
    if (length(lineage[mother_row]) - pos) >= 6
        for id in range(start = pos, stop = length(lineage[mother_row]) - 6, step = 3) # the division for both daughters happens at the same time because in lineage we store the time when the mother cell divides so the next id in lineage is different
            if lineage[mother_row][id] != lineage[mother_row][id + 3]
                time_distance_mother = round(lineage[mother_row][id + 4] - time_div, digits = 3)
                return ((lineage[mother_row][id + 3], time_distance_mother), (lineage[mother_row][id + 3] + 1, time_div)) #(id_daughter_1, time of mother), id_daughter_2)
            end
        end
    else 
        if lineage[mother_row][pos] != lineage[mother_row][pos + 3]
            time_distance_mother = round(lineage[mother_row][pos + 4] - time_div, digits = 3)
            return ((lineage[mother_row][pos + 3], time_distance_mother), (lineage[mother_row][pos + 3] + 1, time_div)) #(id_daughter_1, time of mother), (id_daughter_2, division_time))
        end
    end
    return nothing # the cell root_id did not divide
end


function daughters_id_2(lineage::Vector{Vector{Float64}}, mother_row::Int, pos::Int; mother_time::Vector= [])
    if mother_time == []
        mother_time = [lineage[mother_row][pos + 1]]
    end
    
    time_div = lineage[mother_row][pos + 1]
 # the division for both daughters happens at the same time because in lineage we store the time when the mother cell divides so the next id in lineage is different
    if lineage[mother_row][pos] != lineage[mother_row][pos + 3]
        time_distance_mother = round(lineage[mother_row][pos + 4] - mother_time[1], digits = 3)
        if length(lineage[mother_row]) > 6
            return ((lineage[mother_row][pos + 3], time_distance_mother), (lineage[mother_row][pos + 3] + 1, mother_time[1])) #(id_daughter_1, time of mother), id_daughter_2)
        else 
            return (lineage[mother_row][pos + 3], time_distance_mother)
        end
    end
    return time_div # the cell root_id did not divide it,s a transition not a division
end

function first_div(cell_lineage::Vector, root_id::Float64)
    ind_root = findfirst(x -> x == root_id, cell_lineage)
    for i in range(start = ind_root, stop = length(cell_lineage) - 3, step = 3)
        if cell_lineage[ind_root] != cell_lineage[i + 3]
            return (ind_root, i + 3)
        end
    end
    return nothing
end

function newick_lineage(cells_data::DataFrame, lineage::Vector{Vector{Float64}}, root_id::Float64, newick::String; found_row::Int = 0, start_cell::Bool = false)
    type_dict = Dict()

    if found_row != 0
        mother_row = found_row
    else
        mother_row = find_root(lineage, root_id, cell_count = 1)[1]
    end
    
    ind = first_div(lineage[mother_row], root_id)

    if (length(newick) == 0) & (start_cell == false)
        newick = string(Int(root_id))*";"
        type_dict[string(Int(root_id))] = lineage[mother_row][ind[1] + 2]
        type_dict[string(Int(lineage[mother_row][ind[2]]))] = lineage[mother_row][ind[2] + 2]
        start_pos = ind[2] - 3
        
    elseif (length(newick) == 0) & (start_cell == true)
        newick = "("*string(Int(lineage[mother_row][ind[2]]))*":"*string(round(lineage[mother_row][ind[2] + 1] - lineage[mother_row][ind[1] + 1], digits=3))*")"*string(Int(root_id))*";"
        type_dict[string(Int(root_id))] = lineage[mother_row][ind[1] + 2]
        type_dict[string(Int(lineage[mother_row][ind[2]]))] = lineage[mother_row][ind[2] + 2]
        start_pos = ind[2]
    else
        start_pos = 4
    end

    for ids in range(start = start_pos, stop = length(lineage[mother_row]) - 3, step = 3)
        root_id = string(Int(lineage[mother_row][ids]))
        daughters = daughters_id(lineage, mother_row, ids)
        
        if daughters != nothing
            type_dict[string(Int(daughters[1][1]))] = lineage[mother_row][ids + 5]
            row_daughter_2 = find_root(lineage, daughters[2][1])

            if row_daughter_2 != []
                time_daughter_2 = ":"*string(round(lineage[row_daughter_2[1]][2] - daughters[2][2], digits = 3))
                type_dict[string(Int(daughters[2][1]))] = lineage[row_daughter_2[1]][3]
                
            else
                cell_row_2 = find_cell(cells_data, daughters[2][1])
                time_daughter_2 = ":"*string(round(lineage[cell_row_2][length(lineage[cell_row_2]) - 1] - daughters[2][2], digits = 3))
                type_dict[string(Int(lineage[cell_row_2][length(lineage[cell_row_2]) - 2]))] = lineage[cell_row_2][length(lineage[cell_row_2])]
                
            end
            
            root_in_newick = findfirst(root_id, newick)[1]
            
            newick = newick[1:(root_in_newick - 1)]*"("*string(Int(daughters[1][1]))*":"*string(daughters[1][2])*","*string(Int(daughters[2][1]))*time_daughter_2*")"*newick[root_in_newick:length(newick)]
            if row_daughter_2 != []
                
                newick,sub_type_dict= newick_lineage(cells_data, lineage, daughters[2][1], newick; found_row = row_daughter_2[1], start_cell = false)
                merge(type_dict, sub_type_dict)
            end
        end
    end
    return newick, type_dict
end


function newick_lineage_2(cells_data::DataFrame, lineage::Vector{Vector{Float64}}, root_id::Float64, newick::String; found_row::Int = 0, start_cell::Bool = false)
    type_dict = Dict()

    if found_row != 0
        mother_row = found_row
    else
        mother_row = find_root(lineage, root_id, cell_count = 1)[1]
    end
    
    ind = first_div(lineage[mother_row], root_id)

    if (length(newick) == 0) & (start_cell == false)
        newick = string(Int(root_id))*";"
        type_dict[string(Int(root_id))] = lineage[mother_row][ind[1] + 2]
        type_dict[string(Int(lineage[mother_row][ind[2]]))] = lineage[mother_row][ind[2] + 2]
        start_pos = 1
        if length(lineage[mother_row]) == 6
            time_daughter = string(round(lineage[mother_row][5]- lineage[mother_row][2], digits = 3))
            newick = "("*string(Int(lineage[mother_row][4]))*":"*time_daughter*")"*string(Int(root_id))*";"
            return newick, type_dict
        end
        
    elseif (length(newick) == 0) & (start_cell == true)
        newick = "("*string(Int(lineage[mother_row][ind[2]]))*":"*string(round(lineage[mother_row][ind[2] + 1] - lineage[mother_row][ind[1] + 1], digits=3))*")"*string(Int(root_id))*";"
        type_dict[string(Int(root_id))] = lineage[mother_row][ind[1] + 2]
        type_dict[string(Int(lineage[mother_row][ind[2]]))] = lineage[mother_row][ind[2] + 2]
        start_pos = ind[2]
    else
        start_pos = 1
        
    end

    cell_times = []
    for ids in range(start = start_pos, stop = length(lineage[mother_row]) - 4, step = 3)
        root_id = string(Int(lineage[mother_row][ids]))
        daughters = daughters_id_2(lineage, mother_row, ids; mother_time = cell_times)
        print()
        
        if typeof(daughters) == Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}
            type_dict[string(Int(daughters[1][1]))] = lineage[mother_row][ids + 5]
            row_daughter_2 = find_root(lineage, daughters[2][1])
            

            if row_daughter_2 != []
                time_daughter_2 = ":"*string(round(lineage[row_daughter_2[1]][2] - daughters[2][2], digits = 3))
                type_dict[string(Int(daughters[2][1]))] = lineage[row_daughter_2[1]][3]
                
            else
                cell_row_2 = find_cell(cells_data, daughters[2][1])
                time_daughter_2 = ":"*string(round(lineage[cell_row_2][length(lineage[cell_row_2]) - 1] - daughters[2][2], digits = 3))
                type_dict[string(Int(lineage[cell_row_2][length(lineage[cell_row_2]) - 2]))] = lineage[cell_row_2][length(lineage[cell_row_2])] 

            end
            
            root_in_newick = findfirst(root_id, newick)[1]
            
            newick = newick[1:(root_in_newick - 1)]*"("*string(Int(daughters[1][1]))*":"*string(daughters[1][2])*","*string(Int(daughters[2][1]))*time_daughter_2*")"*newick[root_in_newick:length(newick)]
            if row_daughter_2 != []
                print("tteeessss")
                newick,sub_type_dict= newick_lineage(cells_data, lineage, daughters[2][1], newick; found_row = row_daughter_2[1], start_cell = false)
                merge(type_dict, sub_type_dict)
            end    
        else
            push!(cell_times, daughters)
        end
    end
    return newick, type_dict
end

function stem_newick_lineage(cells_data::DataFrame, lineage::Vector{Vector{Float64}}, root_id::Float64, stem_cells_id)
    if root_id in stem_cells_id
        mother_rows = find_root(lineage, root_id; cell_count = 2)
        newick_1,type_dict_1 = newick_lineage_2(cells_data, lineage, root_id, ""; found_row = mother_rows[1], start_cell = true)[1:2]
        
        newick_2,type_dict_2 = newick_lineage_2(cells_data, lineage, root_id, ""; found_row = mother_rows[2], start_cell = true)[1:2]
        print(newick_2)
        print(type_dict_2)
        root_in_newick_2 = findfirst(")"*string(Int(root_id))*";", newick_2)[1]
        root_in_newick_1 = findfirst(")"*string(Int(root_id))*";", newick_1)[1]

        full_newick = newick_1[1:(root_in_newick_1 - 1)]*","*newick_2[2:(root_in_newick_2 - 1)]*newick_1[root_in_newick_1:length(newick_1)]
        type_dict = merge(type_dict_1, type_dict_2)

        return full_newick,type_dict
    else
        return "not a stem cell"
    end

end

###########
# TO TEST
###########
"""
    collection_from_lineage(
        data : DataFrame containing data of cells at the end of the run, 
        deaths : DataFrame containing data of the cells which died during the run , 
        last_id : last given id
        )

Return adjacent matrix from the data (asymetrix matrix which cooresponds to a directed graph) and type of each cell as a dictionnary (key = cell id and value = type of the cell)

"""

function collection_from_lineage(data::DataFrame, last_id::Int, death_file::String)
    #data is a DataFrame each row corresponds to an alive cell with colunms id, type, lineage
    #last_id is the last ID of all cells which have existed
    #:lineage is in the form of an array recording every transition and division in the form of [id, step, type, id, step, type, ...]
    
    deaths = DataFrame(CSV.File(death_file))
    cell_data = vcat(data, deaths)
    adj_matrix = zeros(Int, last_id, last_id)
    division_type = Dict{Int, Int}(cell_data[!, :id][i] => cell_data[!, :type][i] for i in 1:nrow(cell_data))
    for cell in 1:nrow(cell_data)

        if isa(cell_data[cell, :lineage], String)
            cell_data[cell, :lineage] = parse.(Int, split(chop(cell_data[cell, :lineage], head = 1, tail = 1), ","))
        end
        lineage_data = cell_data[cell, :lineage]
        L = length(lineage_data)

        for i in range(start = 4, stop = L, step = 3)
            if lineage_data[i] != lineage_data[i - 3]
                adj_matrix[lineage_data[i - 3], lineage_data[i]] = 1
                division_type[lineage_data[i - 3]] = lineage_data[i - 1]
            end
        end
    end

    return adj_matrix, division_type
end


#################################################################
# TO REWRITE FOR THE NEW WAY DATA ARE STORED
#################################################################

"""
    get_adj_list_all_cells(
        data : DataFrame containing data of cells at the end of the run,
        last_id : last given id,
        death_file : Path of file containing CSV file containing data of the cells which died during the run 
        )

Return adjacent matrix from the data (asymetrix matrix which cooresponds to a directed graph) and type of each cell as a dictionnary (key = cell id and value = type of the cell)

"""
function get_adj_list_all_cells(data::DataFrame, last_id::Int, death_file::String)
    deaths = DataFrame(CSV.File(death_file))
    cell_data = vcat(data, deaths)
    adj_list = [[] for _ in 1:last_id]
    division_type = Dict{Int, Int}(cell_data[!, :id][i] => cell_data[!, :type][i] for i in 1:nrow(cell_data))
    
    for cell in 1:nrow(cell_data)
        if isa(cell_data[cell, :lineage], String)
            cell_data[cell, :lineage] = parse.(Int, split(chop(cell_data[cell, :lineage], head = 1, tail = 1), ","))
        end
        lineage_data = cell_data[cell, :lineage]
    
        sort!(cell_data,:id)
    
        L_cell = length(lineage_data)
    
        if lineage_data[L_cell - 2] != cell_data[cell, :id]
            push!(adj_list[lineage_data[L_cell - 2]], cell_data[cell, :id])
            division_type[lineage_data[L_cell - 2]] = lineage_data[L_cell]
        end

        for i in range(start = L_cell - 2, stop = 4, step = -3)
            print(i)
            if lineage_data[i] != lineage_data[i - 3]  
                if (lineage_data[i] in adj_list[lineage_data[i - 3]]) == false
                    push!(adj_list[lineage_data[i - 3]], lineage_data[i])
                    division_type[lineage_data[i - 3]] = lineage_data[i - 1]
                end
            end
        end
    end
    return adj_list, division_type
end



function get_adj_list_all_cells2(data, last_id, death_file)
    deaths = DataFrame(CSV.File(death_file))
    cell_data = vcat(data, deaths)
    adj_list = [[] for _ in 1:last_id]
    division_type = Dict{Int, Int}(cell_data[!, :id][i] => cell_data[!, :type][i] for i in 1:nrow(cell_data))
    
    for cell in 1:nrow(cell_data)
        if isa(cell_data[cell, :lineage], String)
            cell_data[cell, :lineage] = parse.(Int, split(chop(cell_data[cell, :lineage], head = 1, tail = 1), ","))
        end
        lineage_data = cell_data[cell, :lineage]
    
        sort!(cell_data,:id)
    
        L_cell = length(lineage_data)
    
        for i in range(start = L_cell - 2, stop = 4, step = -3)
            print("  ", i)
            if lineage_data[i] != lineage_data[i - 3] 
                
                if (lineage_data[i] in adj_list[lineage_data[i - 3]]) == false
                    
                    push!(adj_list[lineage_data[i - 3]], lineage_data[i])
                    division_type[lineage_data[i - 3]] = lineage_data[i - 1]
                else
                    break
                end
            end
        end

        if lineage_data[L_cell - 2] != cell_data[cell, :id]
            push!(adj_list[lineage_data[L_cell - 2]], cell_data[cell, :id])
            division_type[lineage_data[L_cell - 2]] = lineage_data[L_cell]
        end
        
        
    end

    return adj_list, division_type
end



function store_cellcount(file::String, output_file::String, count_func::Function, index::Vector{Int}, last_time_step::Int, nbr_cell::Int; precision::Float64 = 1.0)

    DictofDict = Dict{Int, Dict}(x => Dict() for x in index)
    for i in index

        data_file = file*"/output_run_"*string(i)
        death_file = file*"/output_death_run_"*string(i)
        
        DictofDict[i] = count_func(data_file, death_file, [1, 2, 3, 4, 5], last_time_step, nbr_cell, precision) 
    end
    
    save(output_file,"DictofDict"*string(index[1])*"_"*string(index[length(index)]), DictofDict)
    return (output_file, "DictofDict"*string(index[1])*"_"*string(index[length(index)]))

end