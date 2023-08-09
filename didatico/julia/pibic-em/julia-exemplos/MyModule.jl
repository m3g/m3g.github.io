
module MyModule

  export Test

  # Definition of the structure (type)

  mutable struct Test
    a :: Int64
    b :: Float64
  end

  # Contructor (initialization of the type)

  Test(;a=1,b=2.0) = Test(a,b)
 
end

