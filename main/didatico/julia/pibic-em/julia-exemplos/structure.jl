
#push!(LOAD_PATH,"./")
#using MyModule

#or

include("MyModule.jl")
using .MyModule

x = Test(2,3.0)
println(x)

