#
# Simplex minimization
#
function simplex(f,x0,niter)
  # Initial point: x0 contains 3 vectors of dimension 2
  x = copy(x0)
  # Initial function values
  fx = zeros(3)
  for i in 1:3
    fx[i] = f(x[i])
  end
  xtemp = zeros(2)
  ftemp = 0.
  xav = zeros(2)
  xtrial = zeros(2)
  println(" Initial points: ")
  for i in 1:3
    println(x[i]," ",fx[i])
  end
  # Convergence criterium desired
  convcrit = 1.e-10
  # Main interation
  for iter in 1:niter
    println(" ------ ITERATION: ", iter)
    # Order the points from best to worst
    order = sort([1,2,3],by=i->fx[i])
    x = x[order]
    fx = fx[order] 
    # Check convergence
    if (fx[3]-fx[2] < convcrit) && (fx[3]-fx[1] < convcrit)
      println(" Precision reached. ")
      println(" Best point found: ", x[1], " f = ", fx[1])
      return x[1], fx[1]
    end
    # Compute averge of best points
    @. xav = 0.5*(x[1]+x[2])
    # Compute trial point
    @. xtrial = x[3] + 2*(xav-x[3])
    ftrial = f(xtrial) 
    # If ftrial is better than fx[3], replace point 3 with trial point
    if ftrial < fx[3]
      fx[3] = ftrial
      @. x[3] = xtrial
      println(" Accepted point: ", x[3]," f = ", fx[3])
    else
      println(" Function increased. Trying line search. ")
      # Try up to 10 different points in the 
      # direction x[3]+gamma*(xtrial-x[3])
      for j in 1:10
        @. xtemp = x[3] + rand() * (xtrial - x[3])
        ftemp = f(xtemp)
        if ftemp < fx[3]
          fx[3] = ftemp
          @. x[3] = xtemp
          println("   Line search succeeded at trial ", j)
          println("   New point: ", x[3], " f = ", fx[3])
          break # exits from line search loop
        end
      end
      # If the line search didn't find a better point, stop
      if ftemp > fx[3]
        println(" End of search. ")
        println(" Best point found: ", x[1], " f = ", fx[1])
        return x[1], fx[1]
      end
    end
  end
  println(" Maximum number of trials reached. ")
  println(" Best point found: ", x[1], " f = ", fx[1])
  return x[1], fx[1]
end
