using SyntheticExpressionMixtures
using Documenter

DocMeta.setdocmeta!(SyntheticExpressionMixtures, :DocTestSetup,
                    :(using SyntheticExpressionMixtures); recursive=true)

makedocs(;
         modules=[SyntheticExpressionMixtures],
         authors="Chris Damour",
         sitename="SyntheticExpressionMixtures.jl",
         format=Documenter.HTML(;
                                canonical="https://damourChris.github.io/SyntheticExpressionMixtures.jl",
                                edit_link="main",
                                assets=String[],),
         pages=["Home" => "index.md"],)

deploydocs(;
           repo="github.com/damourChris/SyntheticExpressionMixtures.jl",
           devbranch="main",)
