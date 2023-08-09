
mutable struct Test

  a :: Int64
  b :: Vector{Float64}

end 

# constructors

Test() = Test(0,[0.])

Test(;a=0,b=[0.]) = Test(a,b)

# setting up

x = Test(2,[0.]) # default

y = Test() # defined

z = Test(a=3) # defined



