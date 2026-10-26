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
        "Getting Started" => "getting_started.md",
        "Your Data" => [
            "Survey and model files" => "data/files.md",
            "Building a mesh"        => "data/mesh.md",
            "Control files"          => "data/control_files.md",
        ],
        "Forward Modelling" => "forward.md",
        "Inversion" => [
            "Running an inversion" => "inversion/running.md",
            "Gauss-Newton"         => "inversion/gauss_newton.md",
            "Neural (INR)"         => "inversion/inr.md",
            "Choosing between them" => "inversion/choosing.md",
        ],
        "Plots" => "plots.md",
        "Worked Examples" => "examples.md",
        "API" => "api.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/JuliaGeophysics/GravityGeophysics.jl.git",
    devbranch = "main",
    push_preview = true,
)
