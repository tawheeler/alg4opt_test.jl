if isa(Pkg.installed("Vec"), Void)
    Pkg.clone("https://github.com/tawheeler/Vec.jl.git")
end
if isa(Pkg.installed("ExprRules"), Void)
    Pkg.clone("https://github.com/sisl/ExprRules.jl.git")
end

Pkg.checkout("Weave")
