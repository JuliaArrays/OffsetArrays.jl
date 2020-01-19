using Documenter
using OffsetArrays

makedocs(
    sitename = "OffsetArrays",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = ["index.md", "internals.md", "reference.md"],
    modules = [OffsetArrays],
    doctestfilters = [r"at \./.*", r"at /home.*", r"top-level scope.*", r"\[\d*\]\s*$"],   # for backtraces
)

deploydocs(
    repo = "github.com:JuliaArrays/OffsetArrays.jl.git"
)
