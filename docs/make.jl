using Documenter, DocumenterVitepress
using GravityGeophysics

makedocs(;
    sitename = "GravityGeophysics.jl",
    authors  = "JuliaGeophysics community, Pankaj K Mishra, and contributors",
    format   = DocumenterVitepress.MarkdownVitepress(;
        repo       = "github.com/JuliaGeophysics/GravityGeophysics.jl",
        devbranch  = "main",
        devurl     = "dev",
        # full URL with https://, otherwise the host is taken as part of the base path
        deploy_url = "https://juliageophysics.com/GravityGeophysics.jl",
        description = "A software repository for gravity geophysics research and applications",
        # dev is the only published version, so let search engines index it
        noindex_non_stable = false,
    ),
    # nested lists become the dropdown menus of the top navigation bar
    pages = [
        "Home" => "index.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/JuliaGeophysics/GravityGeophysics.jl.git",
    devbranch = "main",
    push_preview = true,
)
