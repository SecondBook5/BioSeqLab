import random


class GibbsSampler:
    """
    GibbsSampler is an iterative algorithm to find a collection of k-mers that
    represent motifs in a set of DNA sequences. It randomly selects a k-mer
    from each sequence and iteratively updates one k-mer at a time based on a
    profile matrix, aiming to improve the motifs.

    It starts with randomly chosen k-mers in each of t DNA sequences, and in
    each iteration, it randomly selects an integer i between 1 and t and then
    randomly changes only a single k-mer Motif. The process is repeated for N iterations.
    """

    def __init__(self, dna_list, k, t, N, repeats=20):
        """
        Initializes the GibbsSampler object with the given parameters.

        Args:
        - dna_list (list): A list of DNA sequences.
        - k (int): The length of the motifs to be searched for.
        - t (int): The number of DNA sequences in the DNA list.
        - N (int): The number of iterations.
        - repeats (int): The number of times to repeat the algorithm to find the best motifs.
        """
        self.dna_list = dna_list
        self.k = k
        self.t = t
        self.N = N
        self.repeats = repeats

    def gibbs_sampler(self):
        """
        Performs the Gibbs sampling algorithm.

        Returns:
        - best_motifs (list): The best collection of motifs found by the algorithm.
        """
        best_motifs_score = self.k * self.t
        best_motifs = None
        # Repeat the Gibbs sampling algorithm for a specified number of times
        for _ in range(self.repeats):
            bms, bm = self.gibbs_sampler_cycle()
            # Update best motifs if a better set is found
            if bms < best_motifs_score:
                best_motifs = bm
                best_motifs_score = bms
        return best_motifs

    def gibbs_sampler_cycle(self):
        """
        Performs one cycle of the Gibbs sampling algorithm.

        Returns:
        - best_motifs_score (int): The score of the best motifs found in this cycle.
        - best_motifs (list): The best collection of motifs found in this cycle.
        """
        dna_length = len(self.dna_list[0])
        # Initialize motifs with random k-mers from each DNA sequence
        motifs = [self.dna_list[random.randint(
            0, dna_length-self.k):random.randint(0, dna_length-self.k)+self.k] for _ in range(self.t)]
        best_motifs = motifs[:]
        best_motifs_score = self.score(best_motifs)
        # Iterate through N iterations of updating motifs
        for _ in range(self.N):
            i = random.randint(0, self.t-1)
            # Generate profile matrix from motifs excluding the i-th sequence
            profile = self.profile_laplace(
                [motif for index, motif in enumerate(motifs) if index != i])
            # Update the i-th motif based on the profile
            motifs[i] = self.profile_most_probable_kmer(
                self.dna_list[i], profile)
            # Calculate the score of the updated motifs
            motifs_score = self.score(motifs)
            # Update best motifs if the score improves
            if motifs_score < best_motifs_score:
                best_motifs_score = motifs_score
                best_motifs = motifs[:]
        return best_motifs_score, best_motifs

    def score(self, motifs):
        '''
        Calculates the score of the given list of motifs.

        The score is calculated based on the number of nucleotides that do not match the most common nucleotide
        at each position across all motifs.

        Args:
        - motifs (list): A list of motifs (k-mers).

        Returns:
        - score (int): The score of the motifs.
        '''
        # Count occurrences of each nucleotide in each column of the motifs
        counts = [{nucleotide: column.count(
            nucleotide) for nucleotide in 'ACGT'} for column in zip(*motifs)]
        # Calculate the score based on the counts
        return sum(len(motifs) - max(counts[i].values()) for i in range(len(counts)))

    def profile(self, motifs):
        '''
        Constructs the profile matrix of the DNA list motifs.

        The profile matrix represents the probability of each nucleotide at each position in the motifs.

        Args:
        - motifs (list): A list of motifs (k-mers).

        Returns:
        - profile (list of lists): The profile matrix.
        '''
        # Count occurrences of each nucleotide in each column of the motifs and normalize
        counts = [{nucleotide: column.count(
            nucleotide) / len(column) for nucleotide in 'ACGT'} for column in zip(*motifs)]
        return counts

    def profile_most_probable_kmer(self, dna, profile):
        '''
        Finds the profile most probable k-mer for the given input data.

        Args:
        - dna (str): A DNA sequence.
        - profile (list of lists): The profile matrix.

        Returns:
        - most_probable (str): The most probable k-mer in the sequence.
        '''
        max_probability = -1
        most_probable = ''
        # Iterate through all possible k-mers in the given DNA sequence
        for i in range(len(dna) - self.k + 1):
            kmer = dna[i:i+self.k]
            probability = 1
            # Calculate the probability of the k-mer using the profile
            for j, nucleotide in enumerate(kmer):
                probability *= profile[j][nucleotide]
            # Update the most probable k-mer if the probability is higher
            if probability > max_probability:
                max_probability = probability
                most_probable = kmer
        return most_probable

    def profile_laplace(self, motifs):
        '''
        Constructs the Laplace profile matrix of the DNA list motifs.

        Laplace smoothing is applied to avoid zero probabilities.

        Args:
        - motifs (list): A list of motifs (k-mers).

        Returns:
        - profile (list of lists): The Laplace profile matrix.
        '''
        # Count occurrences of each nucleotide in each column of the motifs with Laplace smoothing
        counts = [{nucleotide: (column.count(nucleotide) + 1) / (len(column) + 4)
                   for nucleotide in 'ACGT'} for column in zip(*motifs)]
        return counts
