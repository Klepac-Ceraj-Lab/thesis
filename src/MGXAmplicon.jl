module MGXAmplicon

export
    norm_taxon


const taxlevels = (
    d = :domain,
    k = :kingodom,
    p = :phylum,
    c = :class,
    o = :order,
    f = :family,
    g = :genus,
    s = :species,
    t = :subspecies
)

# Many taxon names are not _actually_ classified to the level suggested.
# eg `Desulfovibrionaceae_unclassified` is listed in genera,
# but only has a family classifier.

unclassified_patterns = (
    r"^Candidatus",
    r"unclassified",
    r"uncultured",
    r"[Ii]ncertae_[Ss]edis",
    r"Clade_I+",
    r"_sp_",
    r"_CAG_",
    r"\d+",
    r"_bacterium$"
)

function norm_taxon(tax::AbstractString; splitter='|', kind=:mgx)
    taxsplit = split(replace(tax, r"[\[\]]"=> ""), splitter)
    if length(taxsplit) == 1 && any(==(lowercase(tax)), ("unassigned", "unclassified", "unknown")) 
        return (taxlevel=:UNKNOWN, taxname=tax)
    end
    kind == :amplicon ? norm_amp_taxon(taxsplit) : norm_taxon(taxsplit)
end

function norm_taxon(taxvec::AbstractVector)
    isempty(taxvec) && throw(ArgumentError("Something went wrong parsing taxon: Vector is empty!"))
    tax = strip(pop!(taxvec))
    
    taxletter = Symbol(first(tax))
    !in(taxletter, keys(taxlevels)) && return norm_taxon(taxvec)
    
    taxlevel = taxlevels[taxletter]
    taxname = replace(tax, r"[kpcofgs]__"=> "")
    return any(p-> occursin(p, taxname), unclassified_patterns) ?
                norm_taxon(taxvec) :
                (; taxlevel, taxname)
end

function norm_amp_taxon(taxvec::AbstractVector)
    taxvec = strip.(taxvec)
    tax = last(taxvec)
    any(p-> occursin(p, replace(tax, r"^[dkpcofgs]__"=>"")), unclassified_patterns) && return norm_amp_taxon(taxvec[1:end-1])
    
    if length(taxvec) >= 7
        sp = taxvec[7]
        # species that end with `sp` or don't start with a capital letter
        # are not acutally classified to the species level
        if any(pat-> occursin(pat, sp), [r"_sp.$", r"^s__[^A-Z]"]) || endswith(sp, "_x")
            # recurse up to genus
            return norm_amp_taxon(taxvec[1:6])
        end
    elseif length(taxvec) == 6
        ge = taxvec[6]
        ge = replace(ge, r"\[|\]"=>"")
        !occursin(r"^g__[A-Z]", ge) && return norm_amp_taxon(taxvec[1:5])
        if occursin(r"_group$", ge)
            m = match(r"g__([A-Z][a-z]+_[a-z]+)_group$", ge)
            if !isnothing(m)
                taxvec[6] = replace(ge, r"_[a-z]+_group$"=>"")
                push!(taxvec, "s__$(m.captures[1])")
                return norm_amp_taxon(taxvec)
            else
                ge = replace(ge, r"_[a-z]+_group$"=>"")
            end
        end
        occursin(r"\d", ge) && return norm_amp_taxon(taxvec[1:5])
        taxvec[6] = ge
    elseif length(taxvec) == 5
        fa = taxvec[5]

    end

    taxlevel = taxlevels[Symbol(first(last(taxvec)))]
    taxname = join(taxvec, ";")
    
    return (; taxlevel, taxname)
end


end # module MGXAmplicon