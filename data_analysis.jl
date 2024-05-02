
##########################################################################################################################################
# Plotting data as graph
##########################################################################################################################################


###########
# TO TEST
###########

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
    
        sort!(cell_data, :id)
    
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
    
        sort!(cell_data, :id)
    
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

function cell_count(data::DataFrame, death_file::String, cell_type::Vector{Int}, last_time_step::Int)
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
                if string(id)*"_"*string(type) in cell_to_count
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
