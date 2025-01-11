from collections import defaultdict
from tqdm import tqdm
from .encoding import encode_kmer, canonical_kmer_bits, generate_all_canonical_kmers_bits

def compute_kmer_signature_bits(sequence, k=2):
    """
    compute the k-mer signature vector for a sequence
    uses bitwise encoding to handle k-mers efficiently
    """
    if len(sequence) < k:
        raise ValueError("Sequence length must be at least k.")

    # dictionary to store frequencies (default value is 0)
    kmer_counts = defaultdict(int)

    # encode the entire sequence into bits
    encoded_sequence = encode_kmer(sequence)

    # mask to extract k-mers
    mask = (1 << (2 * k)) - 1  # mask for the last 2k bits
    for i in range(len(sequence) - k + 1):
        kmer = (encoded_sequence >> (2 * (len(sequence) - k - i))) & mask
        canonical = canonical_kmer_bits(kmer, k)
        kmer_counts[canonical] += 1

    # generate the signature vector
    possible_kmers = generate_all_canonical_kmers_bits(k)
    signature_vector = [kmer_counts[kmer] for kmer in possible_kmers]

    return signature_vector

def signature_for_all_genes(X, k):
    """compute k-mer signatures for all sequences"""
    return [compute_kmer_signature_bits(seq, k) for seq in tqdm(X)]
