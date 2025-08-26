using Documenter, RheoJL

makedocs(
    sitename = "RheoJL Documentation",
    modules = [RheoJL],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
    ],
)