use bio::io::fasta;
use pyo3::class::PyMappingProtocol;
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::PyIterProtocol;
use rand::distributions::{Distribution, Uniform};
use rand::rngs::SmallRng;
use rand::SeedableRng;
use std::fs::File;
use std::io::BufReader;

#[pyclass]
/// Iterator from an Indexed Fasta that samples each record `n_samples` times
/// returning slices of `slice_size`.
/// Length of the sequence has to be strictly greater than `slice_size`
/// to avoid unbalancing since `IndexFastaIterator` is used for training.
pub struct IndexFastaIterator {
    index: fasta::IndexedReader<BufReader<File>>,
    rng: SmallRng,
    /// size of returned slices
    pub slice_size: u64,
    /// number of slices per fasta record
    pub n_samples: usize,
    cycle: bool,
    i: usize,
    len: Option<usize>,
}

impl IndexFastaIterator {
    pub fn new(
        fasta: BufReader<File>,
        fai: fasta::Index,
        slice_size: u64,
        n_samples: usize,
        cycle: bool,
    ) -> Self {
        let mut index_iterator = Self {
            index: fasta::IndexedReader::with_index(fasta, fai),
            rng: SmallRng::from_entropy(),
            slice_size,
            n_samples,
            i: 0,
            len: None,
            cycle,
        };
        // compute records len just once
        index_iterator.records_len();
        index_iterator
    }

    /// Advance to next record.
    fn inner_inc(&mut self) -> Option<()> {
        self.i += 1;
        self.index.fetch_by_rid(self.i - 1, 0, 0).ok()
    }

    // we could avoid inner_inc and get_interval by relying of fetch_by_rid, but
    // we can avoid indexing the same sequence this way
    fn get_interval(&mut self, start: u64, stop: u64) -> Option<Vec<[u8; 4]>> {
        self.index.start = Some(start);
        self.index.stop = Some(stop);
        Some(
            self.index
                .read_iter()
                .ok()?
                .map(|base| match base {
                    Ok(65) => [1, 0, 0, 0], // A
                    Ok(67) => [0, 1, 0, 0], // C
                    Ok(71) => [0, 0, 1, 0], // G
                    Ok(84) => [0, 0, 0, 1], // T
                    _ => [0, 0, 0, 0],
                })
                .collect(),
        )
    }
}

#[pymethods]
impl IndexFastaIterator {
    /// Number of records that satisfy the sequence size minimum
    /// self if mutable because the solution is cached internally
    fn records_len(&mut self) -> usize {
        if let Some(length) = self.len {
            length
        } else {
            let len = self
                .index
                .index
                .sequences()
                .iter()
                .filter(|fasta::Sequence { len: x, .. }| x > &self.slice_size)
                .count()
                * self.n_samples;
            self.len = Some(len);
            len
        }
    }

    /// Get a sequence. Instead of retrieving n_samples, it will retrieve one
    /// sample but taking into account that the length is [`records_len`] * n_samples
    fn get_idx(&mut self, i: usize) -> PyResult<Vec<[u8; 4]>> {
        // usize division is floor division
        let record_index = i / self.__len__();
        // save the fetch operation if we are already in that record
        if (self.i != record_index) & (self.i < self.__len__())
            || (self.index.read(&mut Vec::new()).is_err())
        {
            self.index
                .fetch_by_rid(record_index, 0, 0)
                .map_err(|_| PyIndexError::new_err(format!("{} out of bounds", i)))?;
            self.i = 0;
        }
        if let Some(length) = self.index.inner_len() {
            let between = Uniform::new(0, length - self.slice_size - 1);
            let idx = between.sample(&mut self.rng);
            self.get_interval(idx, idx + self.slice_size)
                .ok_or_else(|| {
                    PyIndexError::new_err(format!(
                        "{} record is malformed; indexed length is wrong.",
                        record_index
                    ))
                })
        } else {
            Err(PyIndexError::new_err(format!(
                "{} record is malformed; it has no defined length",
                record_index
            )))
        }
    }
}

#[pyproto]
impl PyMappingProtocol for IndexFastaIterator {
    /// Actual length of python iterator.
    fn __len__(&self) -> usize {
        self.len.unwrap()
    }
}

#[pyproto]
impl PyIterProtocol for IndexFastaIterator {
    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<Self>) -> Option<Vec<Vec<[u8; 4]>>> {
        slf.next()
    }
}

impl Iterator for IndexFastaIterator {
    type Item = Vec<Vec<[u8; 4]>>;

    fn next(&mut self) -> Option<Self::Item> {
        // get internal next record
        match self.inner_inc() {
            Some(_) => (),
            None => {
                if self.cycle {
                    self.i = 0;
                    return self.next();
                } else {
                    return None;
                }
            }
        }
        match self.index.inner_len() {
            Some(length) if length <= self.slice_size => self.next(),
            None => self.next(),
            Some(length) => {
                // sample N times the current record
                let between = Uniform::new(0, length - self.slice_size - 1);
                (0..self.n_samples)
                    .map(|_| {
                        let idx = between.sample(&mut self.rng);
                        self.get_interval(idx, idx + self.slice_size)
                    })
                    .collect()
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.index.len()))
    }
}
