try Pkg.installed("Vec"); catch
    Pkg.clone("https://github.com/tawheeler/Vec.jl.git")
end
try Pkg.installed("ExprRules"); catch
    Pkg.clone("https://github.com/sisl/ExprRules.jl.git")
end

Pkg.checkout("Weave")
