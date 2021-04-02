---
author:
- Kevin Bonham
- PhD
title: RESONANCE amplicon updated
---

# RESONANCE amplicon updated

## Running the pipeline

Using QIIME.jl

```sh
qiimejl -n test -i v4v5_seqs -o ./ --threads 16 -v
```

### Massage dataframe

```julia
using CSV, DataFrames, Chain


features = CSV.read("/augusta/staff/kevin/echo/echo_amplicon/pipeline/dada2/feature-table.tsv", DataFrame, header=2, delim='\t')
taxa = CSV.read("/augusta/staff/kevin/echo/echo_amplicon/pipeline/taxa/taxonomy.tsv", DataFrame, delim='\t')

rename!(features, "#OTU ID"=>"feature")
rename!(taxa, "Feature ID"=>"feature")

labeled = leftjoin(select(taxa, ["feature", "Taxon"]), features, on="feature")
# some columns don't have any taxonomic assigment, remove them
labeled = labeled[!, [!(eltype(labeled[!, col]) <: Union{Missing, <:Number}) || 
                        sum(labeled[!, col]) > 0 
                        for col in names(labeled)]]

# convert to relative abundance
select!(labeled, [:Taxon => (x -> identity(x)), 
                 [C => c -> c ./ sum(c) for C in names(labeled, r"^[CM]")]...],
                 renamecols=false)
any(isnan, [sum(labeled[:, i]) for i in 2:ncol(labeled)-2])
```

Now, we need to massage the labels a bit,
since we have things like `s__Prevotella_sp.` listed as a species,
and things like ` g__[Ruminococcus]_gnavus_group` not resolved to species.

So we need to consolidate a bit.

```julia
function taxon(levels::AbstractVector)
    if length(levels) >= 7
        sp = levels[7]
        sp = replace(sp, r"\[|\]"=>"")
        # species that end with `sp` or don't start with a capital letter
        # are not acutally classified to the species level
        if any(pat-> occursin(pat, sp), [r"_sp.$", r"^s__[^A-Z]"])
            # recurse up to genus
            return taxon(levels[1:6])
        else
            return (name = replace(sp, "s__"=>""), level = :species)
        end
    elseif length(levels) == 6
        ge = levels[6]
        ge = replace(ge, r"\[|\]"=>"")
        !occursin(r"^g__[A-Z]", ge) && return (name=join(levels, "; "), level=:family_up)
        if occursin(r"_group$", ge)
            m = match(r"g__([A-Z][a-z]+_[a-z]+)_group$", ge)
            if !isnothing(m)
                push!(levels, "s__$(m.captures[1])")
                return taxon(levels)
            else
                ge = replace(ge, r"_group$"=>"")
            end
        end
        return occursin(r"\d", ge) ? (name=join(levels, "; "), level=:family_up) :
                                     (name = split(ge, "_", keepempty=false)[2], level = :genus)
    else
        return (name=join(levels, "; "), level=:family_up)
    end
end


function taxon(longstring::AbstractString)
    levels = strip.(split(longstring, ';'))
    return taxon(levels)
end

labeled = hcat(labeled, DataFrame(taxon.(labeled.Taxon)))

    
```


```julia
labeled.genus = [row.level in (:genus, :species) ? String(split(row.name, "_")[1]) : "UNCLASSIFIED" for row in eachrow(labeled)]
labeled.species = [row.level != :species ? "UNCLASSIFIED" : row.name for row in eachrow(labeled)]

ge = @chain labeled begin
    groupby(:genus)
    combine([n => sum for n in names(labeled, r"^[CM]")], renamecols=false)
end


sp = @chain labeled begin
    groupby(:species)
    combine([n => sum for n in names(labeled, r"^[CM]")], renamecols=false)
end

CSV.write("/home/kevin/repos/danielle-thesis/paper-taxonomic-levels/resonance_dada2_genera.csv", ge)
CSV.write("/home/kevin/repos/danielle-thesis/paper-taxonomic-levels/resonance_dada2_species.csv", sp)

```

### Get some statistics

```julia
using Statistics

sp_totals = sum.(eachcol(filter(:species=> !=("UNCLASSIFIED"), sp)[!, r"^[CM]"]));
mean(sp_totals)
median(sp_totals)

ge_totals = sum.(eachcol(filter(:genus => !=("UNCLASSIFIED"), ge)[!, r"^[CM]"]));
mean(ge_totals)
median(ge_totals)
```
