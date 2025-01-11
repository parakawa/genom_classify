from itertools import product

def encode_kmer(kmer):
    """convert a k-mer to its integer representation"""
    encoding = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}
    result = 0
    for base in kmer:
        result = (result << 2) | encoding[base]
    return result

def decode_kmer(encoded_kmer, k):
    """convert an encoded k-mer back to its string representation"""
    decoding = ['A', 'C', 'G', 'T']
    kmer = []
    for _ in range(k):
        base = encoded_kmer & 0b11  # take the last 2 bits
        kmer.append(decoding[base])
        encoded_kmer >>= 2
    return ''.join(reversed(kmer))

def reverse_complement_bits(kmer, k):
    """compute the reverse complement of a k-mer encoded as an integer"""
    mask = (1 << (2 * k)) - 1  # mask to keep only 2k bits
    complement = 0
    for _ in range(k):
        complement = (complement << 2) | ((kmer & 0b11) ^ 0b11)  # complement the base
        kmer >>= 2
    return complement & mask  # ensure there are no extra bits

def canonical_kmer_bits(kmer, k):
    """return the canonical k-mer in its integer form"""
    rev_comp = reverse_complement_bits(kmer, k)
    return min(kmer, rev_comp)

def generate_all_canonical_kmers_bits(k):
    """generate all canonical k-mers of length k as integers"""
    canonical_kmers = set()
    for kmer_tuple in product([0b00, 0b01, 0b10, 0b11], repeat=k):
        kmer = 0
        for base in kmer_tuple:
            kmer = (kmer << 2) | base
        canonical_kmers.add(canonical_kmer_bits(kmer, k))
    return sorted(canonical_kmers)


if False:
    sequence = "ACGTTCGACG"
    k = 3
    signature = compute_kmer_signature(sequence, k)
    print(f"k-mer signature vector (k={k}): {signature}")
    canonical_kmers = generate_all_canonical_kmers(k)
