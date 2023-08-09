
let Program

  g = -9.8
  x0 = 0.
  y0 = 0.

  dt = 0.001
 
  v0_x = 10.
  v0_y = 10.

  x = [ 0. for i in 1:2000 ]
  y = [ 0. for i in 1:2000 ]
  
  for i in 1:2000
    t = i*dt
    x[i] = x0 + v0_x*t
    y[i] = y0 + v0_y*t + (g/2.)*t^2
  end

  using PyPlot
  plot(x,y)
  savefig("ex1b.png")

end


