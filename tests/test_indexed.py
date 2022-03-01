from lcf import get_fasta_indexed


def test_fasta_iterator_index(genome_fasta: str, genome_index: str):
    for i, seqs in enumerate(
        get_fasta_indexed(genome_fasta, genome_index, 68, 5, False)
    ):
        for seq in seqs:
            assert [1, 0, 0, 0] in seq
            assert [0, 1, 0, 0] in seq
            assert [0, 0, 1, 0] in seq
            assert [0, 0, 0, 1] in seq
    assert i + 1 == 2


def test_fasta_cycle_index(genome_fasta: str, genome_index: str):
    for i, seqs in enumerate(
        get_fasta_indexed(genome_fasta, genome_index, 60, 2, True)
    ):
        print(i)
        if i > 30:
            break

    assert i == 31


def test_fasta_length(genome_fasta: str, genome_index: str):
    assert get_fasta_indexed(genome_fasta, genome_index, 0, 5, False).records_len() == 3
