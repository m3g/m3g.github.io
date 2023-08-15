#!/home/leandro/programs/julia-latest/bin/julia

using HTTP

HTTP.serve() do request::HTTP.Request

  @show request
  HTTP.Response("Hello")

end
