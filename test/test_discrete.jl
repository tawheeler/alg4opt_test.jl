let
    @test is_totally_unimodular(Float64[ 1 -1  0  0 -1;
                                        -1  1 -1  0  0;
                                         0 -1  1 -1  0;
                                         0  0 -1  1 -1;
                                        -1  0  0 -1  1])
    @test is_totally_unimodular(Float64[ 1  1  1  1  1;
                                         1  1  1  0  0;
                                         1  0  1  1  0;
                                         1  0  0  1  1;
                                         1  1  0  0  1])
    @test is_totally_unimodular(Float64[-1 -1  0  0  0  1;
                                         1  0 -1 -1  0  0;
                                         0  1  1  0 -1  0;
                                         0  0  0  1  1 -1])
    @test !is_totally_unimodular(Float64[-1 -1  0  0  0  1;
                                          1  0 -1 -1  0  0;
                                          0  1  1  0 -1  0;
                                          0  0  0  1  2 -1])
    @test !is_totally_unimodular(Float64[-1 -1  0  0  0  1;
                                          1  0  1  1  0  0;
                                          0  1  1 -1 -1  0;
                                          0  0  0  1  1 -1])

    A = Float64[ 1 -1  0  0 -1;
                -1  1 -1  0  0;
                 0 -1  1 -1  0;
                 0  0 -1  1 -1;
                -1  0  0 -1  1]
    b = [1.0, 2, -1, 2, 3]
    c = [2, 1, 3, 4, 5]
    IP = IntegerLinearProgram(A, b, c)
    @test is_totally_unimodular(A, b, c)

    A[1] = 2
    @test !is_totally_unimodular(A, b, c)

    A[1] = 1
    b[1] = 1.5
    @test !is_totally_unimodular(A, b, c)

    A = [0.5 -0.5 1; 2 0.5 -1.5]
    b = [5/2, -3/2]
    c = [2, 1, 3]
    IP = IntegerLinearProgram(A, b, c)
    @test !is_totally_unimodular(A, b, c)

    @test round_ip(IP) == [1,0,2]

    x = cutting_plane(IP)
    @test x == [1,2,3]

    srand(0)
    @test_throws Exception cutting_plane(IntegerLinearProgram(rand(10,10), rand(10), rand(10)))

    x, y = minimize_lp_and_y(IP)
    @test isapprox(x[1], 0.818, atol=1e-3)
    @test isapprox(x[2], 0.000, atol=1e-3)
    @test isapprox(x[3], 2.091, atol=1e-3)
    @test isapprox(y, dot(c, x), atol=1e-10)

    x = branch_and_bound(IP)
    @test x == [1,2,3]

    padovan = [1, 1, 1, 2, 2, 3, 4, 5, 7, 9, 12, 16, 21, 28, 37, 49, 65, 86, 114, 151, 200, 265]
    @test [padovan_topdown(n) for n in 0 : length(padovan)-1] == padovan
    @test [padovan_bottomup(n) for n in 0 : length(padovan)-1] == padovan

    @test knapsack([1], [1], 1) == [1]
    @test knapsack([1], [1], 2) == [1]
    @test knapsack([1], [5], 1) == [0]
    @test knapsack([92,57,49,68,60,43,67,84,87,72],
                   [23,31,29,44,53,38,63,85,89,82], 165) == convert(BitVector,
                   [ 1, 1, 1, 1, 0, 1, 0, 0, 0, 0]) # https://people.sc.fsu.edu/~jburkardt/datasets/knapsack_01/knapsack_01.html

    graph = DiGraph(4)
    add_edge!(graph, 1, 2)
    add_edge!(graph, 2, 3)
    add_edge!(graph, 3, 1)
    add_edge!(graph, 3, 4)
    add_edge!(graph, 2, 4)
    add_edge!(graph, 1, 4)
    add_edge!(graph, 2, 1)
    add_edge!(graph, 3, 2)
    add_edge!(graph, 1, 3)
    add_edge!(graph, 4, 3)
    add_edge!(graph, 4, 2)
    add_edge!(graph, 4, 1)

    srand(0)
    lengths = Dict((1,2)=>4.0, (2,3)=>20, (3,1)=>7, (1,4)=>15, (2,4)=>2, (3,4)=>60,
                   (2,1)=>4.0, (3,2)=>20, (1,3)=>7, (4,1)=>15, (4,2)=>2, (4,3)=>60)

    x = ant_colony_optimization(graph, lengths)
    @test x == [1,3,2,4]
end