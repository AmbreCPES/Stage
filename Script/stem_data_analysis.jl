include("ABMStem_Run.jl")

count_file = "/home/ajaeger/Documents/githubrepos/Stage/Stem_Param_0/Data/count_act_run_"
ind_start = 70
ind_end = 130

full_type0 = []
full_type1 = []
full_type2 = []

first_dif_cell = []
prolif_dif_cell = []

for ind in ind_start:ind_end
    data = readdlm(count_file*string(ind))
    push!(full_type0, data[:, 1])
    push!(full_type1, data[:, 2])
    push!(full_type2, data[:, 3])

    push!(first_dif_cell, findfirst(x -> x > 0, data[:, 3]))
    push!(prolif_dif_cell, linregress(200:500, data[:, 3][200:500]))
end


###  mean plot ###
#p1 = plot(1:length(full_type0[1]),[mean(full_type0), mean(full_type1), mean(full_type1) + mean(full_type0)], ylims= (1,22000), label= ["Quiescent HSCs" "Activated HSCs" "Total HSCs"], title = "Evolution of the number of HSCs in the niche")
#p2 = plot(1:length(full_type0[1]), mean(full_type2), xaxis = "Time (in days)", label= "Differentiated cells" , color = :purple, title = "Evolution of the number of cells exiting the niche")

#plot(p1, p2, layout=(2, 1), legend=:outertopright, titlefont = font(12), xguidefont=font(9), yguidefont=font(9), display = true)

#println(first_dif_cell)

###  Parameter exploration div  ###

p1 = plot(1:length(full_type0[1]),full_type0, color=:RdPu_7, line_z=(70:130)')
p2 = plot(1:length(full_type0[1]), full_type1, color=:RdPu_7, line_z=(70:130)')
p3 = plot(1:length(full_type0[1]), full_type0 .+ full_type1, color=:RdPu_7, line_z=(70:130)')
p4 = plot(1:length(full_type0[1]), full_type2, color=:RdPu_7, line_z=(70:130)', xaxis = "Time (in days)")
p5 = scatter(ind_start:ind_end, [LinearRegression.slope(prolif_dif_cell[x])[1] for x in eachindex(prolif_dif_cell)], color = :magenta, mode = "markers", title = "Differentiation rate in function of activation time", yaxis = "Differentiation rate", xaxis = "Mean activation time (in days)")

explo = plot(p1, p2, p3, p4, p5, layout=(5, 1), size = (700, 1000), titlefont = font(12), xguidefont=font(9), yguidefont=font(9),legend=false, display = true)
#hline!(explo[5], [mean([LinearRegression.slope(prolif_dif_cell[x])[1] for x in eachindex(prolif_dif_cell)])], color=:black)

savefig(explo, "/home/ajaeger/Documents/githubrepos/Stage/Stem_Param_0/Analysis/act_explo")


#for differentiation time of 104 we find approximately the differentiation rate expected in the paper (1/22)
