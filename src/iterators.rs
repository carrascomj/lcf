use bio::io::fasta;
use pyo3::prelude::*;
use pyo3::PyIterProtocol;
use rand::distributions::{Distribution, Uniform};
use rand::rngs::SmallRng;

#[pyclass]
pub struct FastaIterator {
    pub(crate) iter: fasta::Records<std::io::BufReader<std::fs::File>>,
    pub(crate) path: Option<String>,
    pub(crate) rng: SmallRng,
    pub(crate) slice_size: usize,
}

impl Iterator for FastaIterator {
    type Item = (Vec<[u8; 4]>, String);

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.iter.next();
        if let Some(Ok(rec)) = record {
            let seq = rec.seq();
            let length = seq.len();
            // ignore sequences that has a lower length
            println!("Length {}, slice_size {}", &length, &self.slice_size);
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
            Some((slice, rec.desc().unwrap_or_default().to_string()))
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

#[pyproto]
impl PyIterProtocol for FastaIterator {
    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<Self>) -> Option<(Vec<[u8; 4]>, String)> {
        slf.next()
    }
}

#[pyclass]
pub struct ReFastaIterator {
    pub(crate) iter: fasta::Records<std::io::BufReader<std::fs::File>>,
    pub(crate) path: Option<String>,
    pub(crate) rng: SmallRng,
    pub(crate) slice_size: usize,
    pub(crate) n_samples: usize,
}

#[pyproto]
impl PyIterProtocol for ReFastaIterator {
    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<Self>) -> Option<Vec<Vec<[u8; 4]>>> {
        slf.next()
    }
}

impl Iterator for ReFastaIterator {
    type Item = Vec<Vec<[u8; 4]>>;

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.iter.next();
        if let Some(Ok(rec)) = record {
            let seq = rec.seq();
            let length = seq.len();
            let between = match &self.slice_size.cmp(&length) {
                std::cmp::Ordering::Less => Uniform::from(0..length - self.slice_size),
                _ => return self.next(),
            };
            Some(
                (0..self.n_samples)
                    .map(|_| {
                        let index = between.sample(&mut self.rng);
                        return seq[index..(self.slice_size + index - 1)]
                            .iter()
                            .map(|base| match base {
                                65 => [1, 0, 0, 0], // A
                                67 => [0, 1, 0, 0], // C
                                71 => [0, 0, 1, 0], // G
                                84 => [0, 0, 0, 1], // T
                                _ => [0, 0, 0, 0],
                            })
                            .collect::<Vec<[u8; 4]>>();
                    })
                    .collect(),
            )
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
