using Documenter, RheoJL

makedocs(
    sitename = "RheoJL Documentation",
    modules = [RheoJL],
    format = Documenter.HTML(repolink="git@github.com:CVerhoosel/RheoJL.git", edit_link="main"),
    pages = [
        "Home" => "index.md",
    ],
)