function main()

  j :: Int64 = 0
  j = 2^63 - 1
  println(j) 
  println(typemax(j))

  k :: Int64 = j + 1
  println(k)
  
  y :: Float64 = 0.

  y = 1.79769313e308
  println(y)
  y = y + 0.00000001e308
  println(y)

  s = [ 0. , 1. ]

  l :: Int32 = 0
  l = test(l)

  a = 1
  println(typeof(a))
  b = 1.
  println(typeof(b))
  c = 1.e0
  println(typeof(c))
  d = 1.f0
  println(typeof(d))

end

function test(i)
  i :: Int64 = 1
  return i
end





