"""Fixtures for testing the bdino package."""

from os.path import dirname, join

import pytest


@pytest.fixture
def genome_fasta() -> str:
    return join(
        dirname(__file__),
        "data",
        "three_cds.fna",
    )
