function min(x,f,g,precision)
  xbest = copy(x) # Vector to save the best point
  fbest = f(xbest) # Best value up to now
  println(" Initial point: ",xbest," ",fbest)
  xtrial = copy(x)   # Vector of same dimension
  gbest = g(xtrial)  # Compute gradient at initial point
  gnorm = sqrt(gbest[1]^2+gbest[2]^2+gbest[3]^2)
  deltas = 0.1
  while gnorm > precision
    # Move x in the descent direction, with step deltas
    for i in 1:3
      xtrial[i] = xbest[i] - deltas * gbest[i] # Move x in the -f' direction
    end
    # Compute function value at trial point
    ftrial = f(xtrial)
    # If the function decreased, accept trial point and increase step
    if ftrial < fbest
      # Update best point (do NOT use "xbest = xtrial"!)
      for i in 1:3
        xbest[i] = xtrial[i]
      end
      fbest = ftrial
      # Update gradient
      gbest = g(xbest) 
      gnorm = sqrt(gbest[1]^2+gbest[2]^2+gbest[3]^2)
      # Print progress and udpate deltas
      println(" Accepted: ", xbest," ",fbest," ",deltas," ",gbest)
      deltas = deltas * 2
    else
      println(" Not accepted: ", xtrial," ",ftrial," ",deltas,xbest)
      deltas = deltas / 2
    end
  end
  println(" Critical point found. ")
  println(" xbest = ", xbest, " fbest = ", fbest, " g = ", gbest )
  return xbest, fbest
end

f(x :: Vector{Float64}) = x[1]^2 + x[2]^2 + x[3]^2
g(x :: Vector{Float64}) = [ 2*x[1] , 2*x[2] , 2*x[3] ]

x = [ 5.0, 7.0, -3.0 ]
precision = 1.e-8

xmin, fmin = min(x,f,g,precision)
println("xmin = ", xmin," fmin = ", fmin)
