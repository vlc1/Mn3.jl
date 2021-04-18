using Pluto

jl2html(str) = replace(str, r".jl"i => ".html")

function pluto2html(input,
                    output = jl2html(input))
    ss = Pluto.ServerSession()  
    nb = Pluto.SessionActions.open(ss, input; run_async=false)
    html = Pluto.generate_html(nb)
    write(output, html)
end

dir = joinpath(dirname(pwd()), "sample")
files = readdir(dir)
n = length(files)

for file in files
    input = joinpath(dir, file)
    output = joinpath(pwd(), "src", "sample", jl2html(file))

    pluto2html(input, output)
end

open(joinpath("src", "notebook.md"), "a") do io
    for (i, file) in enumerate(files)
        write(io, "1. [Notebook $i](" * joinpath("sample", jl2html(file)) * (i == n ? ").\n\n" : "),\n"))
    end
end

using Documenter, Mn3

makedocs(sitename="Mn3",
         pages = ["Home" => "index.md",
                  "Notebooks" => "notebook.md"])

