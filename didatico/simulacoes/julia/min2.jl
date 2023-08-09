function min2(x0,precision)
  xbest = x0 # Save best point
  fbest = xbest^2 # Best value up to now
  println(" Initial point: ",xbest," ",fbest)
  deltaf = -1.
  deltax = 0.1
  while deltaf < 0. 
    # Move x in the descent direction, with step deltax
    dfdx = 2*xbest # Computing the derivative
    if abs(dfdx) < precision
      println(" Critical point found. ")
      println(" xbest = ", xbest, " fbest = ", fbest, " dfdx = ", dfdx )
      return xbest, fbest
    end
    xtrial = xbest - deltax * dfdx # Move x in the -f' direction
    # Compute function value at trial point
    ftrial = xtrial^2
    # If the function decreased, accept trial point and increase step
    if ftrial < fbest
      xbest = xtrial
      fbest = ftrial
      println(" Accepted: ", xbest," ",fbest," ",deltax," ",dfdx)
      deltax = deltax * 2
    else
      println(" Not accepted: ", xbest," ",fbest," ",deltax," ",dfdx)
      deltax = deltax / 2
    end
  end
end
