function sum_a(a0 :: Float32)
  a = a0 
  for i in 1:10000
    a = a + 1.e-4
  end
  return a
end
