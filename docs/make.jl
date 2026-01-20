import Pkg

# Activate the docs environment (isolated from global)
Pkg.activate(@__DIR__)

try
    # Ensure docs environment uses the local, unregistered package
    Pkg.develop(path=joinpath(@__DIR__, ".."))
    Pkg.develop(url="https://github.com/wangfh5/DataProcessforDQMC.jl.git")
    Pkg.instantiate()
catch e
    @warn "Docs environment setup (Pkg.develop/instantiate) failed" e
end

using Documenter
using BSAHelper

DocMeta.setdocmeta!(BSAHelper, :DocTestSetup, :(using BSAHelper); recursive=true)

makedocs(;
    modules=[BSAHelper],
    authors="ssqc and contributors",
    sitename="BSAHelper.jl",
    format=Documenter.HTML(;
        canonical="https://ssqc.github.io/BSAHelper.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "概览" => "index.md",
        "模块一：bsa_core.jl" => "bsa_core.md",
        "模块二：bsa_bootstrap.jl" => "bsa_bootstrap.md",
    ],
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/ssqc/BSAHelper.jl",
    devbranch="main",
)

# Optional: Start local preview server
# Usage: DOCS_PREVIEW=true julia docs/make.jl
if get(ENV, "DOCS_PREVIEW", "false") == "true"
    build_dir = joinpath(@__DIR__, "build")
    port = parse(Int, get(ENV, "DOCS_PORT", "8275"))

    println("\n" * "="^60)
    println("Starting local documentation preview server...")
    println("URL: http://localhost:$(port)")
    println("Stop: Press Ctrl+C")
    println("="^60 * "\n")

    cd(build_dir) do
        run(`python3 -m http.server $(port)`)
    end
end
