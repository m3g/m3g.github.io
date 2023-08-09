
G = 6.67e-11  # m^3/kg s^-2 
massa_terra = 5.972e24 # em kg
raio_terra = 6.371e6 # em metros

massa_corpo = 1.0 # kg

r = raio_terra + 1.

teta = 30 * ( pi / 180 ) 

x = r * cos(teta)
y = r * sin(teta)

d = sqrt( x^2 + y^2 )

F = G * massa_terra * massa_corpo / d^2 

Fx = -F * cos(teta) # aponta na direção de x negativo
Fy = -F * sin(teta) # aponta na direção de y negativo

# A

println(" F = ", F)
println(" Fx = ", Fx)
println(" Fy = ", Fy)

# B

ax = Fx / massa_corpo
ay = Fy / massa_corpo

println(" ax = ", ax )
println(" ay = ", ay )

# C

x0 = x # calculado acima
y0 = y # calculado acima 
v0x = 0.
v0y = 0.

t = 10 # escolha um tempo

x = x0 + v0x*t + (ax/2)*t^2
y = y0 + v0y*t + (ax/2)*t^2

println(" t = ", t)
println(" x = ", x )
println(" y = ", y )

