def metagenomicsAbundanceNormalizer(otu_abundance_file, 
                                       species_name, 
                                           abundance_criteria):
    """summary_line
    a metagenomics abundance normalizer which will take the abundance 
    OTUs file and gives you a normalized ratio for plotting of the species
    Normally when you plot OTUs, the OTU graph doesnt make sense as there are
    lot of OTUs and the clumpsy graph, this approach, takes a species abundance
    and then divide it by the total number of the OTUs associates with the 
    species and gives you a number. You can plot as many species as you want
    and the analysis will look better and meaningful. 
    
    Keyword arguments:
    argument --
    otu_abundance_file_ : metagenomics abundance OTUs, see the sample file
    species_name_ : species name for which you want to predict the score
    abundance_criteria_ : means filter the OTUs less than those counts.
    """
    
    import pandas as pd
    abundance = int(abundance_criteria)
    taxonomy_read = pd.read_csv(otu_abundance_file, sep = "\t")
    taxonomy_read["split"] = taxonomy_read["Taxonomy"].apply(lambda n: n.split(";"))
    taxonomy_read["taxa_split"] = taxonomy_read["split"].apply(lambda n: list(filter \
                                                    (None,[i.replace("(100)", "").replace("\"","") for i in n])))
    species_filter = taxonomy_read[taxonomy_read["taxa_split"].apply(lambda n: species_name in n)]
    number_filter = species_filter.iloc[::].where(species_filter["Size"] >= abundance).dropna()
    OTU_species_sum = sum(number_filter["Size"])
    OTU_length = len(number_filter["OTU"])
    return OTU_species_sum//OTU_length
