#
# Use the simplex minimizer
#
include("./simplex-struct.jl")

#
# Compute the function value (performs a simulation to do so)
#
function computef(x,data)

  # The constants are given in x 
  k1 = x[1]
  km1 = x[2]

  # Number of data points
  ndata = length(data.Aexp)

  # Vector that will contain the simulation data
  Asim = similar(data.Aexp)

  # Initial concentrations
  Asim[1] = data.CA0
  CB = data.CB0

  # Simulate the reaction using the parameters given
  for i in 2:ndata
    Asim[i] = Asim[i-1] - k1*Asim[i-1]*dt + km1*CB*dt
    CB = CA0 + CB0 - Asim[i]
  end

  # Compute the deviation relative to experimental data
  computef = 0.
  for i in 1:ndata
    computef = computef + ( Aexp[i] - Asim[i] )^2
  end
 
  return computef

end 

# Read experimental data
using DelimitedFiles
data = readdlm("./julia/kineticmodel.dat")
t = data[:,1]
Aexp = data[:,2]
dt = t[2] - t[1] # Time-step
println(" Time-step found: ", dt)

# Initial guess for rate constants, for simplex
# The "experimental" constants in sim2 were k1=0.8 and km1=0.3
k1 = 2.
km1 = 1.
x0 = [ Vector{Float64}(undef,2) for i in 1:3 ]
for i in 1:3
  x0[i][1] = k1 + 0.1*(rand() - 0.5)
  x0[i][2] = km1 + 0.1*(rand() - 0.5)
end

# Initial experimental concentrations
CA0 = 10.
CB0 = 0.

# Define structure
struct Data
  Aexp :: Vector{Float64}
  CA0 :: Float64
  CB0 :: Float64
  dt :: Float64
end
data = Data(Aexp,CA0,CB0,dt)

# Call optimizer
niter = 1000
simplex(computef,x0,niter,data)


