let Sim1
  file = open("sim1.dat","w")
  nsteps = 1000 # Number of steps
  dt = 1.e-1 # Time-step
  k1 = 0.1e0 # Velocity constant
  CA = 10.e0 # Initial concentration
  time = 0.e0
  write(file, "# Time    CA \n")
  for i in 1:nsteps
    CA = CA - k1*CA*dt
    time = time + dt
    write(file, "$time $CA \n")
  end
  close(file)
end
