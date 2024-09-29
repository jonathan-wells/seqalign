use regex::Regex;
use std::io::{Error, ErrorKind};
use std::{collections::HashMap, fs::read_to_string, i32, usize};

#[derive(Debug)]
pub struct ScoringMatrix {
    matrix: HashMap<(char, char), i32>,
}

pub struct AlignmentResult {
    pub s1: Vec<char>,
    pub s2: Vec<char>,
    pub matrix: Vec<Vec<i32>>,
    pub traceback_matrix: Vec<Vec<u8>>,
}

pub struct Aligner {
    scoring_matrix: ScoringMatrix,
    gap_open: i32,
}

impl ScoringMatrix {
    pub fn from_file(filename: &str) -> Result<ScoringMatrix, Error> {
        let data: String = match read_to_string(filename) {
            Ok(data) => data,
            Err(e) => {
                eprintln!("File does not exist: {filename}.");
                return Err(e);
            }
        };
        let lines: Vec<&str> = data.lines().collect();

        // Create HashMap mapping indices in matrix to associated amino acid.
        let aa_map: HashMap<usize, char> =
            lines[1].chars().filter(|c| *c != ' ').enumerate().collect();

        // Create HashMap of amino acid pairs to scores
        let mut matrix: HashMap<(char, char), i32> = HashMap::new();
        let score_re = Regex::new(r"\-*\d+").unwrap();

        let mut i = 0;
        let mut j = 0;
        for line in &lines[2..] {
            let scores: Vec<i32> = score_re
                .find_iter(line)
                .map(|s| s.as_str().parse().unwrap())
                .collect();
            if scores.len() != lines[2..].len() {
                return Err(Error::new(ErrorKind::InvalidData, format!("Dimensions of data are invalid: {filename}")));
            }
            for score in scores {
                matrix.insert((aa_map[&i], aa_map[&j]), score);
                j += 1;
            }
            i += 1;
            j = 0;
        }

        Ok(ScoringMatrix { matrix })
    }

    pub fn get(&self, a: char, b: char) -> i32 {
        self.matrix[&(a, b)]
    }
}

impl AlignmentResult {
    pub fn max_score(&self) -> i32 {
        let idx = self.argmax();
        self.matrix[idx.0][idx.1]
    }

    pub fn alignment(&self) -> (String, String) {
        let mut a1 = String::new();
        let mut a2 = String::new();

        let (mut i, mut j) = self.argmax();
        loop {
            match self.traceback_matrix[i][j] {
                3 => {
                    a1.push('-');
                    a2.push(self.s2[j - 1]);
                    j -= 1
                }
                2 => {
                    a1.push(self.s1[i - 1]);
                    a2.push('-');
                    i -= 1
                }
                1 => {
                    a1.push(self.s1[i - 1]);
                    a2.push(self.s2[j - 1]);
                    i -= 1;
                    j -= 1;
                }
                0 => break,
                _ => panic!(
                    "Encountered unexpected traceback value: {:?}",
                    self.traceback_matrix[i][j]
                ),
            }
        }
        let a1: String = a1.chars().rev().collect();
        let a2: String = a2.chars().rev().collect();
        (a1, a2)
    }

    fn argmax(&self) -> (usize, usize) {
        let mut score: i32 = 0;
        let mut idx: (usize, usize) = (0, 0);

        for i in 1..self.s1.len() + 1 {
            for j in 1..self.s2.len() + 1 {
                if self.matrix[i][j] > score {
                    score += self.matrix[i][j];
                    idx = (i, j);
                }
            }
        }

        idx
    }
}

impl Aligner {
    pub fn new(scoring_matrix: &str, gap_open: i32) -> Result<Self, Error> {
        let scorefile = format!("/Users/jonwells/Projects/fun/seqalign/data/{scoring_matrix}.txt");
        let scoring_matrix = ScoringMatrix::from_file(&scorefile)?;
        let aligner = Aligner {
            scoring_matrix,
            gap_open,
        };
        Ok(aligner)
    }

    pub fn align(&self, s1: String, s2: String) -> AlignmentResult {
        let s1: Vec<char> = s1.chars().collect();
        let s2: Vec<char> = s2.chars().collect();

        let m = s1.len();
        let n = s2.len();

        let mut matrix = vec![vec![0i32; n + 1]; m + 1];
        let mut traceback_matrix = vec![vec![0u8; n + 1]; m + 1];

        for i in 1..m + 1 {
            for j in 1..n + 1 {
                let seqmatch = matrix[i - 1][j - 1] + self.scoring_matrix.get(s1[i - 1], s2[j - 1]);
                let indel_i = matrix[i - 1][j] + self.gap_open;
                let indel_j = matrix[i][j - 1] + self.gap_open;

                let maxval = *(vec![seqmatch, indel_i, indel_j, 0i32]
                    .iter()
                    .max()
                    .expect("Vec will never be empty."));
                matrix[i][j] = maxval;
                traceback_matrix[i][j] = match maxval {
                    x if x == seqmatch => 1,
                    x if x == indel_i => 2,
                    x if x == indel_j => 3,
                    0 => 0,
                    _ => panic!("Unreachable value obtained: {:?}", maxval),
                }
            }
        }

        AlignmentResult {
            s1,
            s2,
            matrix,
            traceback_matrix,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scoring_matrix_symmetry() {
        let blosum62_file = "/Users/jonwells/Projects/fun/seqalign/data/BLOSUM62.txt";
        let scoring_matrix = ScoringMatrix::from_file(blosum62_file).unwrap();
        assert_eq!(scoring_matrix.get('A', 'R'), scoring_matrix.get('R', 'A'));
    }

    #[test]
    fn missing_blosum() {
        let blosum62_file = "/Users/jonwells/Projects/fun/seqalign/data/nofile.txt";
        let sm_err = ScoringMatrix::from_file(blosum62_file).unwrap_err();
        assert_eq!(sm_err.kind(), ErrorKind::NotFound);
    }

    #[test]
    fn bad_blosum() {
        let blosum62_file = "/Users/jonwells/Projects/fun/seqalign/data/BLOSUM62_bad.txt";
        let sm_err = ScoringMatrix::from_file(blosum62_file).unwrap_err();
        assert_eq!(sm_err.kind(), ErrorKind::InvalidData);
    }
}
