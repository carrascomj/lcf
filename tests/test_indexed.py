from lcf import get_fasta_indexed


def test_fasta_iterator_index(genome_fasta: str, genome_index: str):
    i = 0
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
    i = 0
    for i, _ in enumerate(
        get_fasta_indexed(genome_fasta, genome_index, 60, 2, True)
    ):
        print(i)
        if i > 30:
            break

    assert i == 31


def test_fasta_length_is_correct(genome_fasta: str, genome_index: str):
    n_samples = 5
    assert (
        len(get_fasta_indexed(genome_fasta, genome_index, 0, n_samples, False))
        == 3 * n_samples
    )


def test_accessing_is_great(genome_fasta: str, genome_index: str):
    n_samples = 5
    fasta = get_fasta_indexed(genome_fasta, genome_index, 0, n_samples, False)
    assert fasta.get_idx(len(fasta) - 1) == fasta.get_idx(3 * n_samples -1 )
