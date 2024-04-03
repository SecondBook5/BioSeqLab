import random


def score(motifs):
    '''Calculates the score of the given list of motifs.

    The score is calculated as the sum of differences between the most frequent nucleotide 
    in each column of the motif matrix and the other nucleotides.

    Args:
    - motifs (list): A list of motifs (strings).

    Returns:
    - score (int): The score of the motifs.
    '''
    # Create a list of strings, where each string represents a column of the motifs matrix
    columns = [''.join(seq) for seq in zip(*motifs)]
    # Calculate the maximum count of any nucleotide in each column
    max_count = sum([max([c.count(nucleotide)
                    for nucleotide in 'ACGT']) for c in columns])
    # Calculate the score based on the sum of differences between the maximum count and the other counts
    return len(motifs[0]) * len(motifs) - max_count


def profile(motifs):
    '''Constructs the profile of the DNA list motifs.

    The profile is a matrix where each column represents the probability of each nucleotide
    at each position in the motifs.

    Args:
    - motifs (list): A list of motifs (strings).

    Returns:
    - profile (list): The profile matrix.
    '''
    # Create a list of strings, where each string represents a column of the motifs matrix
    columns = [''.join(seq) for seq in zip(*motifs)]
    # Calculate the probability of each nucleotide at each position in the motifs
    return [[float(col.count(nuc)) / float(len(col)) for nuc in 'ACGT'] for col in columns]


def profile_most_probable_kmer(dna, k, profile):
    '''Finds the profile most probable k-mer for the given input data.

    Args:
    - dna (str): A DNA sequence.
    - k (int): The length of the k-mer.
    - profile (list): The profile matrix.

    Returns:
    - most_probable (str): The most probable k-mer in the sequence.
    '''
    # A dictionary relating nucleotides to their position within the profile.
    nuc_loc = {nucleotide: index for index, nucleotide in enumerate('ACGT')}

    # Initialize the maximum probability.
    max_probability = -1

    # Compute the probability of each k-mer, store it if it's currently a maximum.
    for i in range(len(dna)-k+1):
        # Get the current probability.
        current_probability = 1
        for j, nucleotide in enumerate(dna[i:i+k]):
            # Multiply the probability of each nucleotide at each position
            current_probability *= profile[j][nuc_loc[nucleotide]]

        # Check for a maximum.
        if current_probability > max_probability:
            max_probability = current_probability
            most_probable = dna[i:i+k]

    return most_probable


def profile_laplace(motifs):
    '''Constructs the Laplace profile of the DNA list motifs.

    The Laplace profile is a matrix where each column represents the probability of each nucleotide
    at each position in the motifs with Laplace smoothing applied.

    Args:
    - motifs (list): A list of motifs (strings).

    Returns:
    - profile (list): The profile matrix.
    '''
    # Create a list of strings, where each string represents a column of the motifs matrix
    columns = [''.join(seq) for seq in zip(*motifs)]
    # Calculate the Laplace-smoothed probability of each nucleotide at each position in the motifs
    return [[float(col.count(nuc) + 1) / float(len(col) + 4) for nuc in 'ACGT'] for col in columns]


def _gibbs_sampler_cycle(dna_list, k, t, N):
    '''Performs a single cycle of the Gibbs Sampler algorithm.

    Args:
    - dna_list (list): A list of DNA sequences.
    - k (int): The length of the motifs to search for.
    - t (int): The number of sequences in dna_list.
    - N (int): The number of iterations to perform.

    Returns:
    - best_motifs_score (int): The score of the best motifs found.
    - best_motifs (list): The best motifs found.
    '''
    # Calculate the length of the DNA sequences
    dna_length = len(dna_list[0])
    # Generate random start indices for the motifs in each sequence
    random_integers = [random.randint(0, dna_length-k) for i in range(t)]
    # Initialize motifs with random k-mers from each sequence
    motifs = [dna_list[i][r:r+k] for i, r in enumerate(random_integers)]
    # Initialize the best motifs and score with the initial motifs
    best_motifs = list(motifs)
    best_motifs_score = score(best_motifs)
    # Perform N iterations of the Gibbs Sampler algorithm
    for j in range(N):
        # Randomly choose a sequence index to update its motif
        i = random.randint(0, t-1)
        # Construct a profile matrix from the other motifs
        profile = profile_laplace(
            [motif for index, motif in enumerate(motifs) if index != i])
        # Update the motif for sequence i based on the profile
        motifs[i] = profile_most_probable_kmer(dna_list[i], k, profile)
        # Calculate the score of the updated motifs
        motifs_score = score(motifs)
        # Update the best motifs and score if the updated motifs have a lower score
        if motifs_score < best_motifs_score:
            best_motifs_score = motifs_score
            best_motifs = motifs
    return best_motifs_score, best_motifs


def gibbs_sampler(dna_list, k, t, N):
    '''Runs the Gibbs Sampler algorithm to find the best motifs.

    Args:
    - dna_list (list): A list of DNA sequences.
    - k (int): The length of the motifs to search for.
    - t (int): The number of sequences in dna_list.
    - N (int): The number of iterations to perform.

    Returns:
    - best_motifs (list): The best motifs found.
    '''
    # Initialize the best motifs score with the maximum possible score
    best_motifs_score = k * t
    best_motifs = None
    # Repeat the Gibbs Sampler algorithm multiple times to increase the chance of finding the best motifs
    for repeat in range(30):
        # Perform a single cycle of the Gibbs Sampler algorithm
        bms, bm = _gibbs_sampler_cycle(dna_list, k, t, N)
        # Update the best motifs and score if the current cycle produced a better result
        if bms < best_motifs_score:
            best_motifs = bm
            best_motifs_score = bms
    return best_motifs
