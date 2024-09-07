using SynthethicExpressionMixtures
using Documenter

DocMeta.setdocmeta!(SynthethicExpressionMixtures, :DocTestSetup, :(using SynthethicExpressionMixtures); recursive=true)

makedocs(;
    modules=[SynthethicExpressionMixtures],
    authors="Chris Damour",
    sitename="SynthethicExpressionMixtures.jl",
    format=Documenter.HTML(;
        canonical="https://damourChris.github.io/SynthethicExpressionMixtures.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/damourChris/SynthethicExpressionMixtures.jl",
    devbranch="main",
)
