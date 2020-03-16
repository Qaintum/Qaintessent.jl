my_tests = ["test_circuit.jl",
            "test_gates.jl"]

println("Running tests:")
for my_test in my_tests
    include(my_test)
end

