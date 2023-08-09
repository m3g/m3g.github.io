g = -9.8
t = 0.5
x0 = 0.
y0 = 0.

angulo = 30 * ( pi / 180 ) # variar este valor
vmodulo = sqrt( 10^2 + 10^2 )
v0x = vmodulo * cos(angulo)
v0y = vmodulo * sin(angulo)

npassos = 2500
x = [ 0. for i in 1:npassos ]
y = [ 0. for i in 1:npassos ]
for i in 1:npassos
  dt = 0.001
  t = i*dt
  x[i] = x0 + v0x*t
  y[i] = y0 + v0y*t + (g/2)*t^2
  if y[i] <= 0 
    println(" t fim = ", t, " x fim = ", x[i])
    break
  end
end
