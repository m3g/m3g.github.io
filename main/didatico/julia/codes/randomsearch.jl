#
# function that computes the function value
#
function computef(x)
  f = x[1]^2 + x[2]^2
  return f
end
#
# Program randomsearch
#
let RandomSearch
  # Test 10000 points
  x = [ 0. , 0. ]
  xbest = x
  ntrial = 10000
  fbest = computef(x)
  for i in 1:ntrial
    x[1] = -10.e0 + 20.e0*rand(Float64)
    x[2] = -10.e0 + 20.e0*rand(Float64)
    f = computef(x)
    if f < fbest 
      fbest = f
      xbest = x
      println( i, " New best point: ", x, " f = ", f )
    end
  end
  println(" Best point found: ", xbest, " f = ", fbest)
end
