try Pkg.installed("Vec"); catch
    Pkg.clone("https://github.com/tawheeler/Vec.jl.git")
end

Pkg.checkout("Weave")
