using PyPlot

#
# Reading a file using our custom read_data function (see data.jl)
#

push!(LOAD_PATH,"./")
using Data

ndata, x, y = read_data("teste.dat", cols=[1,2], comment="#")

plot(x,y,label="Dados 1")
xlabel("x")
ylabel("y")
legend(loc="lower right")
annotate("Veja aqui",xy=(5,40))
savefig("plot.png")


