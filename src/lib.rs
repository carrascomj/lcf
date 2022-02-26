use bio::io::fasta;
use pyo3::prelude::*;
use pyo3::{wrap_pyfunction, PyIterProtocol};
use rand::distributions::{Distribution, Uniform};
use rand::rngs::SmallRng;
use rand::SeedableRng;

#[pyclass]
pub struct FastaIterator {
    iter: fasta::Records<std::io::BufReader<std::fs::File>>,
    path: Option<String>,
    rng: SmallRng,
    slice_size: usize,
}

#[pyproto]
impl PyIterProtocol for FastaIterator {
    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<Self>) -> Option<Vec<[u8; 4]>> {
        slf.next()
    }
}

/// Retrieve a random slice of each sequence of size slice_size. Additionally,
/// the slice is one hot encoded.
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
            Ok(rec) if rec.seq().len() > slice_size => Some(1),
            _ => None,
        })
        .count()
}

impl Iterator for FastaIterator {
    type Item = Vec<[u8; 4]>;

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.iter.next();
        if let Some(Ok(rec)) = record {
            let seq = rec.seq();
            let length = seq.len();
            // ignore sequences that has a lower length
            let index = match &self.slice_size.cmp(&length) {
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => return self.next(),
                std::cmp::Ordering::Less => {
                    let between = Uniform::from(0..length - self.slice_size);
                    between.sample(&mut self.rng)
                }
            };
            let slice = seq[index..(self.slice_size + index - 1)]
                .iter()
                .map(|base| match base {
                    65 => [1, 0, 0, 0], // A
                    67 => [0, 1, 0, 0], // C
                    71 => [0, 0, 1, 0], // G
                    84 => [0, 0, 0, 1], // T
                    _ => [0, 0, 0, 0],
                })
                .collect();
            Some(slice)
        } else if record.is_some() {
            self.next()
        } else if let Some(data_path) = &self.path {
            // We have reached the end, cycle again
            let reader = fasta::Reader::from_file(data_path).unwrap();
            self.iter = reader.records();
            self.next()
        } else {
            None
        }
    }
}

#[pymodule]
fn lcf(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<FastaIterator>()?;
    m.add_function(wrap_pyfunction!(get_fasta_iterator, m)?)?;
    m.add_function(wrap_pyfunction!(len_fasta_valid, m)?)?;
    Ok(())
}
