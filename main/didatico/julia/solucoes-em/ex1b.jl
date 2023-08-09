g = -9.8
t = 0.5
x0 = 0.
y0 = 0.
v0x = 10.
v0y = 10.
npassos = 2500
x = [ 0. for i in 1:npassos ]
y = [ 0. for i in 1:npassos ]
for i in 1:npassos
  dt = 0.001
  t = i*dt
  x[i] = x0 + v0x*t
  y[i] = y0 + v0y*t + (g/2)*t^2
end
using PyPlot
plot(x,y)
savefig("ex1b.png")
