# Load required packages
using JuliaWebAPI

# Define functions testfn1 and testfn2 that we shall expose
function testfn1(arg1, arg2; narg1="1", narg2="2")
    return (parse(Int, arg1) * parse(Int, narg1)) + (parse(Int, arg2) * parse(Int, narg2))
end

testfn2(arg1, arg2; narg1="1", narg2="2") = testfn1(arg1, arg2; narg1=narg1, narg2=narg2)

# Expose testfn1 and testfn2 via a ZMQ listener
process(
    JuliaWebAPI.create_responder([
        (testfn1, true),
        (testfn2, false)
    ], "tcp://127.0.0.1:9999", true, "")
)
