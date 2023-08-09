
using PyPlot

function main()

  k1 = 1.
  km1 = 2.

  A0 = 1.
  B0 = 0.

  tmax = 2.
  t = 0.
  dt = 0.01

  nsteps = Int64(tmax/dt)
  A = [ 0. for i in 1:nsteps ]
  B = [ 0. for i in 1:nsteps ]
  t = [ 0. for i in 1:nsteps ]
  A[1] = A0
  B[1] = B0
  t[1] = 0.
  for i in 2:nsteps
    t[i] = t[i-1] + dt
    A[i] = A[i-1] + (-A[i-1]*k1 + B[i-1]*km1)*dt
    B[i] = B[i-1] + (A[i-1]*k1 - B[i-1]*km1)*dt
  end

  plot(t,A)
  plot(t,B)
  savefig("cinetica.png")

end 
main()
