let
    grammar = @grammar begin
        R = R + R
        R = R - R
        R = R * R
        R = R / R
        R = log(R) | sin(R) | cos(R) | sqrt(R) | abs(R) | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
    end

    f(node) = begin
        try
            return (eval(node, grammar) - π)^2 + length(node)*1e-5
        catch
            return Inf
        end
    end

    n = 0
    root = RuleNode(13)
    stack = NodeLoc[]
    x_best = RuleNode(0)
    y_best = Inf
    _nested_monte_carlo_expression_discovery(f, grammar, root, stack, n, x_best, y_best)

    n = 0
    root = RuleNode(1, [RuleNode(11), RuleNode(12)])
    stack = NodeLoc[]
    x_best = RuleNode(13)
    y_best = f(x_best)
    _nested_monte_carlo_expression_discovery(f, grammar, root, stack, n, x_best, y_best)

    n = 0
    root = RuleNode(1, [RuleNode(11), RuleNode(0)])
    stack = [NodeLoc(root, 2)]
    x_best = RuleNode(0)
    y_best = Inf
    _nested_monte_carlo_expression_discovery(f, grammar, root, stack, n, x_best, y_best)

    n = 1
    root = RuleNode(1, [RuleNode(11), RuleNode(0)])
    stack = [NodeLoc(root, 2)]
    x_best = RuleNode(0)
    y_best = Inf
    _nested_monte_carlo_expression_discovery(f, grammar, root, stack, n, x_best, y_best)

    n = 3
    root = RuleNode(1, [RuleNode(11), RuleNode(0)])
    stack = [NodeLoc(root, 2)]
    x_best = RuleNode(0)
    y_best = Inf
    _nested_monte_carlo_expression_discovery(f, grammar, root, stack, n, x_best, y_best)

    n = 3
    root = RuleNode(1, [RuleNode(0), RuleNode(0)])
    stack = [NodeLoc(root, 2), NodeLoc(root, 1)]
    x_best = RuleNode(0)
    y_best = Inf
    _nested_monte_carlo_expression_discovery(f, grammar, root, stack, n, x_best, y_best)

    n = 3
    root = RuleNode(1, [RuleNode(0), RuleNode(0)])
    stack = [NodeLoc(root, 2), NodeLoc(root, 1)]
    x_best = RuleNode(0)
    y_best = Inf
    _nested_monte_carlo_expression_discovery(f, grammar, root, stack, n, x_best, y_best)

    nested_monte_carlo_expression_discovery(f, grammar, :R, 1)
    nested_monte_carlo_expression_discovery(f, grammar, :R, 3)

    decisions = [1,2]
    root = complete_expression(decisions, grammar, :R)
    @test root.ind == 1
    @test root.children[1].ind == 2
    eval(root)

    srand(0)
    arr = mutate(IntegerGaussianMutation(1.0), [1,2,3,4])
    @test isa(arr, Vector{Int})
    @test length(arr) == 4

    srand(0)
    @test mutate(GeneDuplication(), [1]) == [1,1]
    mutate(GeneDuplication(), [1,2,3,4,5,6,7,8])

    srand(0)
    arr = mutate(GenePruning(1.0, grammar, :R), [10,1,1,1,1,1,1,1]) == [10]
end
let
    grammar = @grammar begin
        G10 =  R+A # 1
        G25 =  R+A # 2
        G30 =  R+A # 3
        G50 =  R+A # 4
        G60 =  R+A # 5
        G100 = R+A # 6
        A = A+A  # 7
        A = G10  # 8
        A = G25  # 9
        A = G30  # 10
        A = G50  # 11
        A = G60  # 12
        A = G100 # 13
        A = H    # 14
        R = R+R  # 15
        R = G10  # 16
        R = G25  # 17
        R = G30  # 18
        R = G50  # 19
        R = G60  # 20
        R = G100 # 21
        R = E    # 22
        G10 = E
        G25 = E
        G30 = E
        G50 = E
        G60 = E
        G100 = E
        A = E
    end

    radius(sym::Symbol) = parse(Float64, string(sym)[2:end])

    function get_hand_periods!(
        node::RuleNode,
        r::Float64=25.0, # radius of the parent
        t::Float64=10.0, # turn period of the parent
        link::Symbol=:A,
        periods::Vector{Float64}=Float64[]
        )

        typ = return_type(grammar, node)
        rule = grammar.rules[node.ind]
        if rule == :H
            push!(periods, t) # a hand has same turn period as parent
        elseif rule == :E
            # do nothing
        elseif typ == :A || typ == :R
            for c in node.children
                get_hand_periods!(c, r, t, typ, periods) # same properties
            end
        elseif ismatch(r"G\d+", string(typ))
            r2 = radius(typ)
            t2 = link == :R ? -t*r/r2 : t
            for c in node.children
                get_hand_periods!(c, r2, t2, :A, periods)
            end
        end
        return periods
    end

    f  = (node) -> begin
        periods = get_hand_periods!(node)
        if isempty(periods)
            return Inf
        end
        return minimum((1-t)^2 for t in periods) +
               minimum((60-t)^2 for t in periods) +
               minimum((3600-t)^2 for t in periods) +
               length(node)*1e-2
    end

    node = RuleNode(22)
    @test get_hand_periods!(node) == Float64[]

    node = RuleNode(2, [RuleNode(22), RuleNode(8, [RuleNode(1, [RuleNode(21, [RuleNode(6, [RuleNode(21, [RuleNode(6, [RuleNode(22), RuleNode(14)])]), RuleNode(29)])]), RuleNode(29)])])])
    @test get_hand_periods!(node) == [1.0]

    node = RuleNode(2, [RuleNode(22), RuleNode(7, [
        RuleNode(8, [RuleNode(1, [RuleNode(21, [RuleNode(6, [RuleNode(21, [RuleNode(6, [RuleNode(22), RuleNode(14)])]), RuleNode(29)])]), RuleNode(29)])]),
        RuleNode(12, [RuleNode(5, [RuleNode(16, [RuleNode(1, [RuleNode(16, [RuleNode(1, [RuleNode(22), RuleNode(14)])]), RuleNode(29)])]), RuleNode(29)])]),
        ])])

    srand(0)
    max_depth = 10
    population = [rand(RuleNode, grammar, :G25, max_depth) for i in 1 : 100]
    k_max = 100
    S = TournamentSelection(5)
    C = TreeCrossover(grammar, 20)
    M = MultiMutate([TreeMutation(grammar, 0.05),
                     TreePermutation(grammar, 0.05)])

    x = genetic_algorithm(f, population, k_max, S, C, M)
    @test f(x) < 1000.0
end
let
    grammar = @grammar begin
        S = NP * " " * VP
        NP = ADJ * " " * NP
        NP = ADJ * " " * N
        VP = V  * " " * ADV
        ADJ = |(["a", "the", "big", "little", "blue", "red"])
        N = |(["mouse", "cat", "dog", "pony"])
        V = |(["ran", "sat", "slept", "ate"])
        ADV = |(["quietly", "quickly", "soundly", "happily"])
    end

    @test eval(decode([2,10,19,0,6], grammar, :S).node, grammar) == "little dog ate quickly"


    grammar = @grammar begin
        R  =  D * De * P * E
        De  =  D * De | ""
        P  =  "." * D * De | ""
        E  =  "E" * S * D * De | ""
        S  =  "+"|"-"|""
        D  =  "0"|"1"|"2"|"3"|"4"|"5"|"6"|"7"|"8"|"9"
    end

    x = [205, 52, 4, 27, 10, 59, 6]
    str = eval(decode(x, grammar, :R).node, grammar)
    @test str == "4E+8"
end
let
    grammar = @grammar begin
        R = R + R
        R = R * R
        R = F
        R = I
        F = 1.5
        I = 1
        I = 2
        I = error("never choose me")
    end

    probgram = ProbabilisticGrammar(grammar, Dict(:R => [1,1,5,5], :F => [1], :I => [1,1,0]))
    node = RuleNode(1, [RuleNode(3, [RuleNode(5)]), RuleNode(4, [RuleNode(7)])]) # 1.5 + 2
    @test isapprox(probability(probgram, node), (1/12)*(5/12)*(5/12)*(1/2))

    update!(probgram, [node])
    @test probgram.ws[:R] == [1,0,1,1]
    @test probgram.ws[:F] == [1]
    @test probgram.ws[:I] == [0,1,0]

    ppt = PPTNode(grammar)
    @test isempty(ppt.children)
    c = get_child(ppt, grammar, 1)
    @test length(ppt.children) == 1

    node = rand(ppt, grammar, :R)

    @test isa(probability(ppt, grammar, node), Float64)
    @test isa(p_target(ppt, grammar, node, 1.0, 0.5, 0.5, 0.5), Float64)
    @test isa(update!(ppt, grammar, node, 1.0, 0.5, 0.5, 0.5, 0.5), PPTNode)
    @test isa(mutate!(ppt, grammar, node, 0.1, 0.5), PPTNode)
    @test isa(prune!(ppt, grammar), PPTNode)
end