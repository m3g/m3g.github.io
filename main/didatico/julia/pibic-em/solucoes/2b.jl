
G = 6.67e-11  # m^3/kg s^-2 
massa_terra = 5.972e24 # em kg
raio_terra = 6.371e6 # em metros

massa_corpo = 1.0 # kg

r = raio_terra + 1.

teta = 45 * ( pi / 180 ) 

x = r * cos(teta)
y = r * sin(teta)

d = sqrt( x^2 + y^2 )

F = G * massa_terra * massa_corpo / d^2 

println(" F = ", F)
