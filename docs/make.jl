using TimeSeriesTools
using Documenter

DocMeta.setdocmeta!(TimeSeriesTools, :DocTestSetup, :(using TimeSeriesTools); recursive=true)

makedocs(;
    modules=[TimeSeriesTools],
    authors="John Waczak <john.louis.waczak@gmail.com>",
    repo="https://github.com/john-waczak/TimeSeriesTools.jl/blob/{commit}{path}#{line}",
    sitename="TimeSeriesTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://john-waczak.github.io/TimeSeriesTools.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/john-waczak/TimeSeriesTools.jl",
    devbranch="main",
)
