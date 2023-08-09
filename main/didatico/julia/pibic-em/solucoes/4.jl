using PyPlot

let TerraLua

  G = 19.968
  
  raio_terra = 1.
  distancia_lua = 60.
  
  massa_terra = 1.
  massa_lua = 0.0123
  
  npassos = 30 * 24 * 60 # um mes, um passo por minuto

  # Iniciando os vetores de posição, velocidade e tempo

  x = [ 0. for i in 1:npassos ]
  y = [ 0. for i in 1:npassos ] 
  vx = [ 0. for i in 1:npassos ]
  vy = [ 0. for i in 1:npassos ] 
  t = [ 0. for i in 1:npassos ] 

  # posicao e velocidade iniciais

  x[1] = 0.
  y[1] = 60. # raios da terra
 
  vx[1] = 0.565 # raios da terra por hora
  vy[1] = 0.

  t[1] = 0.

  dt = 1. / 60 # um minuto, em horas
  for i in 2:npassos

    t[i] = i*dt
  
    r = sqrt( x[i-1]^2 + y[i-1]^2 )
    F = G * massa_terra * massa_lua / r^2

    Fx = F * ( x[i-1] / r )
    Fy = F * ( y[i-1] / r )

    ax = -Fx / massa_lua
    ay = -Fy / massa_lua

    vx[i] = vx[i-1] + ax*dt
    vy[i] = vy[i-1] + ay*dt
  
    x[i] = x[i-1] + vx[i-1]*dt + (ax/2)*dt^2
    y[i] = y[i-1] + vy[i-1]*dt + (ay/2)*dt^2

    # Para criar o arquivo para o PyMol
    if i % 100 == 0
      println("2")
      println("Lua")
      println("C ",0.,0.,0.)
      println("H ",x[i]," ",y[i]," ",0.)
    end

  end

  plot(x,y)
  savefig("terra_e_lua.png")

end

  

