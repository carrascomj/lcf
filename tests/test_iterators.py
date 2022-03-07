from lcf import get_fasta_iterator, len_fasta_valid


def test_fasta_iterator(genome_fasta: str):
    for i, (seq, desc) in enumerate(get_fasta_iterator(genome_fasta, 65, False)):
        assert [1, 0, 0, 0] in seq
        assert [0, 1, 0, 0] in seq
        assert [0, 0, 1, 0] in seq
        assert [0, 0, 0, 1] in seq
        assert desc
    assert i == 2


def test_fasta_cycle(genome_fasta: str):
    for i, (seq, desc) in enumerate(get_fasta_iterator(genome_fasta, 65, True)):
        assert [1, 0, 0, 0] in seq
        assert [0, 1, 0, 0] in seq
        assert [0, 0, 1, 0] in seq
        assert [0, 0, 0, 1] in seq
        if i > 30:
            break

    assert i == 31


def test_fasta_length(genome_fasta: str):
    assert len_fasta_valid(genome_fasta, 10) == 3
    assert len_fasta_valid(genome_fasta, 90) == 2
    assert len_fasta_valid(genome_fasta, 10000000) == 0
