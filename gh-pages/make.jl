using Pluto

replaceext(str; in = "jl", out = "html") =
    replace(str, Regex(".$(in)", "i") => ".$(out)")

function pluto2html(infile, outfile)
    ss = Pluto.ServerSession()
    nb = Pluto.SessionActions.open(ss, infile; run_async=false)
    html = Pluto.generate_html(nb)
    write(outfile, html)
end

ghpages = pwd()
build = joinpath(ghpages, "build")
root = dirname(ghpages)

innb = joinpath(root, "notebook")
outnb = joinpath(build, "notebook")

mkdir(build)
mkdir(outnb)

files = readdir(innb)

for file in files
    infile = joinpath(innb, file)
    outfile = joinpath(outnb, replaceext(file))

    pluto2html(infile, outfile)
end

open(joinpath(outnb, "index.md"), "w") do io
    write(io, "A list of notebooks:\n\n")

    for (i, file) in enumerate(files)
        write(io, "1. [Notebook $i]("
              * joinpath("sample", replaceext(file))
              * (i == length(files) ? ").\n\n" : "),\n"))
    end
end

