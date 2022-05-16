use bio::io::fasta;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use rand::rngs::SmallRng;
use rand::SeedableRng;

mod indexed_iterators;
mod iterators;

pub use indexed_iterators::IndexFastaIterator;
pub use iterators::{FastaIterator, ReFastaIterator};

/// Retrieve a random slice of each sequence of size `slice_size`. Additionally,
/// the slice is one hot encoded.
///
/// Only sequences greater or equal than `slice_size` will be iterated. This
/// is because the non-indexed iterators are expected to be used for evaluation
/// where sequences should not be discarded.
#[pyfunction]
pub fn get_fasta_iterator(path: &str, slice_size: usize, cycle: bool) -> FastaIterator {
    // TODO: handle this unwrap
    let reader = fasta::Reader::from_file_with_capacity(1000000000, path).unwrap();
    let rng = SmallRng::from_entropy();
    let path = if cycle { Some(path.to_string()) } else { None };
    FastaIterator {
        iter: reader.records(),
        path,
        rng,
        slice_size,
    }
}

/// Retrieve lenght of fasta file in term of valid sequences
#[pyfunction]
pub fn len_fasta_valid(path: &str, slice_size: usize) -> usize {
    // TODO: handle this unwrap
    let reader = fasta::Reader::from_file_with_capacity(1000000000, path).unwrap();
    reader
        .records()
        .filter_map(|rrec| match rrec {
            Ok(rec) if rec.seq().len() >= slice_size => Some(1),
            _ => None,
        })
        .count()
}

/// Retrieve a random slice of each sequence of size slice_size. Additionally,
/// the slice is one hot encoded.
///
/// Only sequences strictly greater than `slice_size` will be iterated to avoid
/// unbalancing of `slice_size` length sequences (which would be repeated
/// `n_samples` times otherwise).
#[pyfunction]
pub fn get_fasta_reiterator(
    path: &str,
    slice_size: usize,
    cycle: bool,
    n_samples: usize,
) -> ReFastaIterator {
    // TODO: handle this unwrap
    let reader = fasta::Reader::from_file_with_capacity(1000000000, path).unwrap();
    let rng = SmallRng::from_entropy();
    let path = if cycle { Some(path.to_string()) } else { None };
    ReFastaIterator {
        iter: reader.records(),
        path,
        rng,
        slice_size,
        n_samples,
    }
}

#[pyfunction]
/// Retrieve an interator that takes `n_samples` of size `n_samples for each
/// record in an indexed FASTA and FASTA.fai files.
///
/// Only sequences strictly greater than `slice_size` will be iterated to avoid
/// unbalancing of strictly equal length sequences (which would be repeated
/// `n_samples` times otherwise).
pub fn get_fasta_indexed(
    path: &str,
    index_path: &str,
    slice_size: u64,
    n_samples: usize,
    cycle: bool,
) -> IndexFastaIterator {
    // TODO: handle this unwrap
    let index = fasta::Index::from_file(&std::path::PathBuf::from(index_path)).unwrap();
    IndexFastaIterator::new(
        std::io::BufReader::new(std::fs::File::open(path).unwrap()),
        index,
        slice_size,
        n_samples,
        cycle,
    )
}

#[pymodule]
fn lcf(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<FastaIterator>()?;
    m.add_class::<ReFastaIterator>()?;
    m.add_function(wrap_pyfunction!(get_fasta_iterator, m)?)?;
    m.add_function(wrap_pyfunction!(get_fasta_reiterator, m)?)?;
    m.add_function(wrap_pyfunction!(len_fasta_valid, m)?)?;
    m.add_function(wrap_pyfunction!(get_fasta_indexed, m)?)?;
    Ok(())
}
