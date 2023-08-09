using PyPlot

let TerraLua

  G = 19.968
  
  raio_terra = 1.
  distancia_lua = 60.
  
  massa_terra = 1.
  massa_lua = 0.0123
  
  npassos = 30 * 24 * 60 # um mes, um passo por minuto

  # Iniciando os vetores de posição, velocidade e tempo

  x_lua = [ 0. for i in 1:npassos ]
  y_lua = [ 0. for i in 1:npassos ] 
  vx_lua = [ 0. for i in 1:npassos ]
  vy_lua = [ 0. for i in 1:npassos ] 

  x_terra = [ 0. for i in 1:npassos ]
  y_terra = [ 0. for i in 1:npassos ] 
  vx_terra = [ 0. for i in 1:npassos ]
  vy_terra = [ 0. for i in 1:npassos ] 

  # tempo
  t = [ 0. for i in 1:npassos ] 

  # posicao e velocidade iniciais

  x_lua[1] = 0.
  y_lua[1] = 60. # raios da terra
 
  vx_lua[1] = 0.565 # raios da terra por hora
  vy_lua[1] = 0.

  x_terra[1] = 0.
  y_terra[1] = 0.
  # para que o centro de massa esteja parado
  vx_terra[1] = -vx_lua[1] * massa_lua / massa_terra
  vy_terra[1] = 0.

  t[1] = 0.

  dt = 1. / 60 # um minuto, em horas
  for i in 2:npassos

    t[i] = i*dt
  
    r = sqrt( ( x_lua[i-1] - x_terra[i-1])^2 + (y_lua[i-1]-y_terra[i-1])^2 )
    F = G * massa_terra * massa_lua / r^2

    Fx = F * ( (x_lua[i-1]-x_terra[i-1]) / r )
    Fy = F * ( (y_lua[i-1]-y_terra[i-1]) / r )

    ax_lua = -Fx / massa_lua
    ay_lua = -Fy / massa_lua

    vx_lua[i] = vx_lua[i-1] + ax_lua*dt
    vy_lua[i] = vy_lua[i-1] + ay_lua*dt
  
    x_lua[i] = x_lua[i-1] + vx_lua[i-1]*dt + (ax_lua/2)*dt^2
    y_lua[i] = y_lua[i-1] + vy_lua[i-1]*dt + (ay_lua/2)*dt^2

    ax_terra = Fx / massa_terra
    ay_terra = Fy / massa_terra

    vx_terra[i] = vx_terra[i-1] + ax_terra*dt
    vy_terra[i] = vy_terra[i-1] + ay_terra*dt
  
    x_terra[i] = x_terra[i-1] + vx_terra[i-1]*dt + (ax_terra/2)*dt^2
    y_terra[i] = y_terra[i-1] + vy_terra[i-1]*dt + (ay_terra/2)*dt^2

    # Para criar o arquivo para o PyMol
    #if i % 100 == 0
    #  println("2")
    #  println("Lua e Terra")
    #  println("H ",x_lua[i]," ",y_lua[i]," ",0.)
    #  println("C ",x_terra[i]," ",y_terra[i]," ",0.)
    #end

  end

  plot(x_lua,y_lua)
  plot(x_terra,y_terra)

  savefig("terra_e_lua.png")

end

  

