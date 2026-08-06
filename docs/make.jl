using Crafts, Roly, Documenter

DocMeta.setdocmeta!(Crafts, :DocTestSetup, :(using Crafts, Roly); recursive=true)

makedocs(; sitename="Crafts.jl",
         modules=[Crafts],
         authors="Maximilian HUEBL <maximilian.huebl@ist.ac.at>",
         checkdocs=:exports,
         format=Documenter.HTML(; canonical="https://goodrichgroup.github.io/Crafts.jl",
                                prettyurls=get(ENV, "CI", nothing) == "true"),
         pages=["Home" => "index.md",
                "Assembly systems" => "assemblysystems.md",
                "Bond potentials" => "potentials.md",
                "Equilibrium" => "equilibrium.md",
                "Kinetics" => "kinetics.md",
                "Stability" => "stability.md",
                "Diffusion and rates" => "rates.md",
                "API reference" => "reference.md"])

deploydocs(; repo="github.com/goodrichgroup/Crafts.jl.git", devbranch="main")
