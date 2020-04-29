samples = readlines("samples.txt");

profiles = String[];

for (root, dirs, files) in walkdir("/lovelace/echo/analysis/engaging/")
    occursin("metaphlan", root) || continue
    filter!(f-> endswith(f, "profile.tsv"), files)
    append!(profiles, joinpath.(root, files))
end

filter!(profiles) do f
    m = match(r"main\/([\w\-]+)_S\d{1,2}_profile\.tsv", f)
    isnothing(m) && error("no match: $f")
    s = String(m.captures[1])
    in(s, samples) || in(replace(s, "-"=>"_"), samples)
end

open("metaphlan_paths.txt", "w") do io
    for p in profiles
        write(io, p * '\n')
    end
end