# Function to be minimized
f(x::Vector{Float64}) = x[1]^2 + x[2]^2
# Minimizer by random search
function randomsearch2(f,ntrial,x0 :: Vector{Float64})
  x = copy(x0)
  xbest = copy(x0)
  fbest = f(xbest)
  for i in 1:ntrial
    x[1] = xbest[1] + 1.e-3*(-1. + 2. * rand())
    x[2] = xbest[2] + 1.e-3*(-1. + 2. * rand())
    fx = f(x)
    if fx < fbest
      fbest = fx
      xbest[1] = x[1]
      xbest[2] = x[2]
      println(i," New best point: ", x," f(x) = ", fx)
    end
  end
  println(" Best point found: ",xbest," f = ", fbest)
end
