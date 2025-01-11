def split_sequence_alternate_overlap(sequence, chunk_size=500, overlap=100):
    fragments = []
    i = 0
    toggle_overlap = True  
    while i <= len(sequence) - chunk_size:
        fragments.append(sequence[i:i + chunk_size])
        if toggle_overlap:
            i += chunk_size - overlap  # overlapping
        else:
            i += chunk_size  # no overlapping
        toggle_overlap = not toggle_overlap 
    return fragments

def split_gens_over_all_genomes(df, chunk_size=500, overlap=100):
    """split all genomes in the dataframe into smaller fragments"""
    X = []
    y = []
    for _, row in df.iterrows():
        sequence = row["genome sequence"]
        organism_name = row["organism name"]

        gens_for_sequence = split_sequence_alternate_overlap(sequence, chunk_size, overlap=overlap)
        X.extend(gens_for_sequence)
        y.extend([organism_name] * len(gens_for_sequence))
    return X, y
