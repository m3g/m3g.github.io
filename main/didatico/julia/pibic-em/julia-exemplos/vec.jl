
function main()

  # explicit type
  x = Vector{Float64}(undef, 2)
  x[1] = 1.
  x[2] = 2.
  println(typeof(x))
  println(x)

  # default (float64)
  y = [ 0. for i in 1:2 ]
  y[1] = 1.
  y[2] = 2.
  println(typeof(y))
  println(y)

  # using float32
  z = [ 0.f0 for i in 1:2 ]
  z[1] = 1.
  z[2] = 2.
  println(typeof(z))
  println(z)

  # Matrix
  A = Matrix{Float64}(undef,2,2)
  for i in 1:2
    for j in 1:2
       A[i,j] = i*j
    end
  end
  println(typeof(A))
  println(A)
  B = A
  C = A * B
  println(C)

  x1 = A*x 
  println(x1)

end ; main()
