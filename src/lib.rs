//! Pairwise Sequence Alignment.
//!
//! Provides a library of tools that can be used to perform pairwise sequence alignment of protein
//! sequences. Currently only implements Smith-Waterman local alignment with affine gap penalty
//! (Gotoh), but this may be extended in the future.

#![allow(dead_code, unused_variables)]
#![feature(portable_simd)]

use regex::Regex;
use std::{
    array,
    io::{Error, ErrorKind},
    collections::HashMap,
    str::FromStr,
    cmp::{max_by_key, min_by_key},
    fs::read_to_string,
    simd::i16x8
};

/// This module provides static strings for each BLOSUM matrix.
pub mod constants;
use crate::constants::*;


/// Holds the scoring system and methods required to generate alignments.
pub struct RognesAligner {
    scoring_matrix: ScoringMatrix,
    gap_open: i16,
    gap_extend: i16,
    gap_open_simd: i16x8,
    gap_extend_simd: i16x8
}

impl RognesAligner {

    /// Creates new `Aligner` instance from given scoring system.
    pub fn new(scoring_matrix: &str, gap_open: i16, gap_extend: i16) -> Result<Self, Error> {
        let scoring_matrix = ScoringMatrix::new(scoring_matrix)?;
        let aligner = RognesAligner {
            scoring_matrix,
            gap_open,
            gap_extend,
            gap_open_simd: i16x8::splat(gap_open),
            gap_extend_simd: i16x8::splat(gap_open),
        };
        Ok(aligner)
    }
    

    /// Returns SIMD-accelerated Smith-Waterman local alignment of two sequences.
    pub fn align(&self, query: &str, target: &str) -> AlignmentResult {
        
        /// Calculates the max of a non-SIMD vector of values required to populate a cell.
        fn calc_max(values: &[i16]) -> i16 {
            let maxval: &i16 = values
                .iter()
                .max()
                .expect("vector will never be empty.");
            *maxval
        }

        let tmp_query: Vec<char> = query.chars().collect();
        let tmp_target: Vec<char> = target.chars().collect();
        
        let query = max_by_key(&tmp_query, &tmp_target, |x| x.len()).to_vec();
        let target = min_by_key(&tmp_query, &tmp_target, |x| x.len()).to_vec();


        // if self.query_profile_flag == false { self.build_query_profile(&query) }
        let m = query.len();
        let n = target.len();
        println!("{m}, {n}");

        let mut h = vec![vec![0i16; n + 1]; m + 1];
        // let mut e = vec![vec![0i16; n + 1]; m + 1];
        // let mut f = vec![vec![0i16; n + 1]; m + 1];

        let mut traceback_matrix = vec![vec![0u8; n + 1]; m + 1];

        // let mut sv_max = i16x8::splat(0);
        
        for k in (1..n-LANES).step_by(LANES) {
            // Upper left triangle
            for j in k..k+LANES {
                for i in 1..LANES+k-j {
                    let seqmatch = h[i - 1][j - 1] + self.scoring_matrix.get(query[i - 1], target[j - 1]);
                    let open_i = h[i - 1][j] - self.gap_open;
                    let open_j = h[i][j - 1] - self.gap_open;
                    // h[i][j] = calc_max(&[seqmatch, open_i, open_j, 0i16]);
                    h[i][j] = 1;
                }
            }
            
            // SIMD block
            for i in LANES..m {
                let h_antidiagonal = i16x8::from_array(array::from_fn(|x| h[(i - 1) - x][(k - 1) + x]));
                let seqmatches = i16x8::from_array(array::from_fn(|x| self.scoring_matrix.get(query[(i - 1) - x], target[(k - 1) + x])));
                let seqmatch_antidiagonal = h_antidiagonal + seqmatches;
                let open_i_antidiagonal = i16x8::from_array(array::from_fn(|x| h[(i - 1) - x][k + x])) - self.gap_open_simd;
                let open_j_antidiagonal = i16x8::from_array(array::from_fn(|x| h[i - x][(k - 1) + x])) - self.gap_open_simd;
                
                let input_matrix = seqmatch_antidiagonal.to_array();
                for x in 0..LANES {
                    h[i-x][k+x] = input_matrix[x];
                }

            }

            // Lower right triangle
            for j in k+1..k+LANES {
                for i in (m-(j-k)..m).rev() {
                    let seqmatch = h[i - 1][j - 1] + self.scoring_matrix.get(query[i - 1], target[j - 1]);
                    let open_i = h[i - 1][j] - self.gap_open;
                    let open_j = h[i][j - 1] - self.gap_open;
                    // h[i][j] = calc_max(&[seqmatch, open_i, open_j, 0i16]);
                    h[i][j] = 1;
                }
            }
        }
        
        // Fill in trailing edges right of SIMD blocks
        for i in 1..m {
            for j in n + 1 -(n % LANES)..n + 1 {
                let seqmatch = h[i - 1][j - 1] + self.scoring_matrix.get(query[i - 1], target[j - 1]);
                let open_i = h[i - 1][j] - self.gap_open;
                let open_j = h[i][j - 1] - self.gap_open;
                // h[i][j] = calc_max(&[seqmatch, open_i, open_j, 0i16]);
                h[i][j] = 1;
            }
        }
        
        for i in 0..m {
            println!("{:?}", h[i]);
        }
        // simd_block(k);
        // lower_right(k, m, match_score, gap_open, &mut h);

        AlignmentResult {
            query,
            target,
            matrix: h,
            traceback_matrix
        }
    }

    fn build_query_profile(&self, query: &str) -> Vec<Vec<i16>> {
        let query: Vec<char> = query.chars().collect();
        let m = query.len();
        let seglen = m + LANES - (m % LANES);
        

        let query_profile = HashMap::<&'static char, i16x8>::new();
        for aa in AA_ALPHABET.iter() {
            let h = 0;
            for i in 0..seglen {
                let j = i;
                let aa_profile = [0i16; LANES];
                for k in 0..LANES {
                    todo!();
                }
            }
        }
        todo!();
    }
}


/// Defines scoring system for vanilla local pairwise alignment.
pub struct Aligner {
    scoring_matrix: ScoringMatrix,
    gap_open: i16,
    gap_extend: i16,
}

impl Aligner {
    /// Creates new `Aligner` instance from given scoring system.
    pub fn new(scoring_matrix: &str, gap_open: i16, gap_extend: i16) -> Result<Self, Error> {
        let scoring_matrix = ScoringMatrix::new(scoring_matrix)?;
        let aligner = Aligner {
            scoring_matrix,
            gap_open,
            gap_extend,
        };
        Ok(aligner)
    }
 
    /// Returns Smith-Waterman local alignment of two sequences.
    pub fn align(&self, query: &str, target: &str) -> AlignmentResult {
        
        /// Calculates the max of a vector of values required to populate a cell.
        fn calc_max(values: &[i16]) -> i16 {
            let maxval: &i16 = values
                .iter()
                .max()
                .expect("vector will never be empty.");
            *maxval
        }
        
        let query: Vec<char> = query.chars().collect();
        let target: Vec<char> = target.chars().collect();

        let m = query.len();
        let n = target.len();

        let mut f = vec![vec![0i16; n + 1]; m + 1];
        let mut g = vec![vec![0i16; n + 1]; m + 1];
        let mut h = vec![vec![0i16; n + 1]; m + 1];

        let mut traceback_matrix = vec![vec![0u8; n + 1]; m + 1];

        for i in 1..=m {
            for j in 1..=n {
                let seqmatch = f[i - 1][j - 1] + self.scoring_matrix.get(query[i - 1], target[j - 1]);
                let open_i = f[i - 1][j] - self.gap_open;
                let open_j = f[i][j - 1] - self.gap_open;
                let extend_i = g[i - 1][j] - self.gap_extend;
                let extend_j = h[i][j - 1] - self.gap_extend;

                g[i][j] = calc_max(&[open_i, extend_i]);
                h[i][j] = calc_max(&[open_j, extend_j]);
                f[i][j] = calc_max(&[seqmatch, g[i][j], h[i][j], 0i16]);
                
                traceback_matrix[i][j] = match f[i][j] {
                    x if x == seqmatch => 1,
                    x if x == g[i][j] => 2,
                    x if x == h[i][j] => 3,
                    0 => 0,
                    _ => panic!("Unexpected value encountered in traceback: {:?}", f[i][j]),
                }
            }
        }

        AlignmentResult {
            query,
            target,
            matrix: f,
            traceback_matrix,
        }
    }
}


/// Record for single fasta sequence
pub struct FastaRecord {
    pub name: String,
    pub seq: String,
}

impl FastaRecord {
    /// Returns length of sequence in record
    pub fn len(&self) -> u32 {
        self.seq.len() as u32
    }
}

impl FromStr for FastaRecord {
    type Err = Error; 
    /// Convert from input fasta string to `FastaRecord`
    fn from_str(s: &str) -> Result<Self, Error> {
        let record_re = Regex::new(r"(?ms)>(\S+)[^\n]*\n(.+)").unwrap();
        let Some(caps) = record_re.captures(s) else {
            return Err(Error::new(ErrorKind::InvalidData, format!("Invalid fasta string.\n{s}")))
        };
        let name = caps[1].to_string();
        let mut seq = caps[2].to_string();
        seq.retain(|c| c != '\n');
        Ok(FastaRecord{ name, seq })
    }
}

/// Given a filename, returns a vector of `FastaRecords`.
pub fn parse_fastafile(filename: &str) -> Result<Vec<FastaRecord>, Error> {
    let data = read_to_string(filename);
    let Ok(data) = data else { 
        return Err(Error::new(ErrorKind::NotFound, format!("{filename} does not exist"))) 
    };
    let fasta_re = Regex::new(r"(?ms)(>[^>]+)").unwrap();
    let mut fastas: Vec<FastaRecord> = Vec::new();
    for cap in fasta_re.captures_iter(&data) {
        let fasta_record = FastaRecord::from_str(&cap[0]);
        match fasta_record {
            Ok(fr) => fastas.push(fr),
            Err(e) => return Err(e)
        };
    }
    Ok(fastas)
}

/// Stores matrices of amino acid pair transition scores (1/2 Bit units).
#[derive(Debug)]
pub struct ScoringMatrix {
    matrix: HashMap<(char, char), i16>,
}

impl ScoringMatrix {

    /// Creates a new scoring matrix from set of available options.
    ///
    /// # Example:
    /// ```
    /// use seqalign::ScoringMatrix;
    /// let scoring_matrix = ScoringMatrix::new("BLOSUM62");
    /// ```
    pub fn new(scoring_matrix: &str) -> Result<Self, Error> {
        match scoring_matrix {
            "BLOSUM30" => ScoringMatrix::from_str(BLOSUM30),
            "BLOSUM35" => ScoringMatrix::from_str(BLOSUM35),
            "BLOSUM40" => ScoringMatrix::from_str(BLOSUM40),
            "BLOSUM45" => ScoringMatrix::from_str(BLOSUM45),
            "BLOSUM50" => ScoringMatrix::from_str(BLOSUM50),
            "BLOSUM55" => ScoringMatrix::from_str(BLOSUM55),
            "BLOSUM60" => ScoringMatrix::from_str(BLOSUM60),
            "BLOSUM62" => ScoringMatrix::from_str(BLOSUM62),
            "BLOSUM65" => ScoringMatrix::from_str(BLOSUM65),
            "BLOSUM70" => ScoringMatrix::from_str(BLOSUM70),
            "BLOSUM75" => ScoringMatrix::from_str(BLOSUM75),
            "BLOSUM80" => ScoringMatrix::from_str(BLOSUM80),
            "BLOSUM85" => ScoringMatrix::from_str(BLOSUM85),
            "BLOSUM90" => ScoringMatrix::from_str(BLOSUM90),
            "BLOSUM100" => ScoringMatrix::from_str(BLOSUM100),
            _ => Err(Error::new(ErrorKind::InvalidInput, format!("{scoring_matrix} not implemented."))),
        }
    }
    
    /// Creates new `ScoringMatrix` instance from a &str.
    fn from_str(input_string: &str) -> Result<Self, Error> {
        let lines: Vec<&str> = input_string.lines().filter(|l| !l.contains('#')).collect(); 
        
        // Create HashMap mapping indices in matrix to associated amino acid.
        let aa_map: HashMap<usize, char> =
            lines[0].chars().filter(|c| *c != ' ').enumerate().collect();

        // Create HashMap of amino acid pairs to scores
        let mut matrix: HashMap<(char, char), i16> = HashMap::new();
        let score_re = Regex::new(r"\-*\d+").unwrap();

        let mut j = 0;
        for (i, line) in lines[1..].iter().enumerate() {
            let scores: Vec<i16> = score_re
                .find_iter(line)
                .map(|s| s.as_str().parse().unwrap())
                .collect();
            if scores.len() != lines[1..].len() {
                return Err(Error::new(ErrorKind::InvalidData, "Dimensions of data are invalid"));
            }
            for score in scores {
                matrix.insert((aa_map[&i], aa_map[&j]), score);
                j += 1;
            }
            j = 0;
        }

        Ok(ScoringMatrix { matrix })
    }
    
    /// Parses scoring matrix file and returns instance of ScoringMatrix.
    ///
    /// # Example:
    /// ```
    /// use seqalign::ScoringMatrix;
    /// let scoring_matrix = ScoringMatrix::from_file("./PAM50.txt");
    /// ```
    pub fn from_file(filename: &str) -> Result<Self, Error> {
        let data: String = match read_to_string(filename) {
            Ok(data) => data,
            Err(e) => {
                eprintln!("File does not exist: {filename}.");
                return Err(e);
            }
        };
        Self::from_str(data.as_str())
    }

    /// Extracts score for given pair of amino acids.
    /// 
    /// The success of this operation is only guaranteed by earlier checks on appropriate input
    /// sequences. 
    pub fn get(&self, a: char, b: char) -> i16 {
        self.matrix[&(a, b)]
    }
}

/// Stores the results of a pairwise sequence alingment.
pub struct AlignmentResult {
    query: Vec<char>,
    target: Vec<char>,
    matrix: Vec<Vec<i16>>,
    traceback_matrix: Vec<Vec<u8>>,
}

impl AlignmentResult {
    /// Returns the score of the optimal alignment.
    pub fn max_score(&self) -> i16 {
        let idx = self.argmax();
        self.matrix[idx.0][idx.1]
    }
    
    /// Returns the pair of sequences corresponding to the optimal alignment.
    ///
    /// This method uses the output and traceback matrices associated with a new `AlignmentResult`
    /// instance to reconstruct the optimal alignment of query and target sequences.
    pub fn alignment(&self) -> (String, String) {
        let mut a1: Vec<char> = Vec::new();
        let mut a2: Vec<char> = Vec::new();

        let (mut i, mut j) = self.argmax();
        loop {
            match self.traceback_matrix[i][j] {
                3 => {
                    a1.push('-');
                    a2.push(self.target[j - 1]);
                    j -= 1;
                }
                2 => {
                    a1.push(self.query[i - 1]);
                    a2.push('-');
                    i -= 1;
                }
                1 => {
                    a1.push(self.query[i - 1]);
                    a2.push(self.target[j - 1]);
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
        let a1: String = a1.iter().rev().collect();
        let a2: String = a2.iter().rev().collect();
        (a1, a2)
    }
    
    /// Return the indices of the top score in the "f" result matrix.
    fn argmax(&self) -> (usize, usize) {
        let mut score: i16 = 0;
        let mut idx: (usize, usize) = (0, 0);

        for i in 1..=self.query.len() {
            for j in 1..=self.target.len() {
                if self.matrix[i][j] >= score {
                    score = self.matrix[i][j];
                    idx = (i, j);
                }
            }
        }

        idx
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scoring_matrix_symmetry() {
        let blosum62_file = "./resources/BLOSUM62";
        let scoring_matrix = ScoringMatrix::from_file(blosum62_file).unwrap();
        assert_eq!(scoring_matrix.get('A', 'R'), scoring_matrix.get('R', 'A'));
    }

    #[test]
    fn missing_blosum() {
        let blosum62_file = "./resources/nofile.txt";
        let sm_err = ScoringMatrix::from_file(blosum62_file).unwrap_err();
        assert_eq!(sm_err.kind(), ErrorKind::NotFound);
    }

    #[test]
    fn bad_blosum() {
        let blosum62_file = "./resources/BLOSUM62_malformatted.txt";
        let sm_err = ScoringMatrix::from_file(blosum62_file).unwrap_err();
        assert_eq!(sm_err.kind(), ErrorKind::InvalidData);
    }

    #[test]
    fn max_score() {
        let aligner = Aligner::new("BLOSUM62", 0, 0).unwrap();
        let alignment = aligner.align(&String::from("RRRRR"), &String::from("RRRRR"));
        assert_eq!(alignment.max_score(), 25);
    }
    
    #[test]
    fn fasta_from_string() {
        let fastastring = ">NP_085137.1 zinc finger protein 436 isoform 1 [Homo sapiens]\nACTG\nGTCA";
        let fasta_record = FastaRecord::from_str(fastastring).unwrap();
        assert_eq!(String::from("NP_085137.1"), fasta_record.name);
        assert_eq!(String::from("ACTGGTCA"), fasta_record.seq);
    }
    
    #[test]
    fn fasta_seqlen() {
        let fastastring = ">record\nACTG\nGTCA";
        let fasta_record = FastaRecord::from_str(fastastring).unwrap();
        assert_eq!(fasta_record.len(), 8);
    }
    
    #[test]
    fn parse_fasta_file() {
        let fasta_records = parse_fastafile("./resources/test_fasta.fa").unwrap();
        assert_eq!(String::from("record1"), fasta_records[0].name);
        assert_eq!(String::from("ACTGGCTA"), fasta_records[0].seq);
        assert_eq!(String::from("record2"), fasta_records[1].name);
        assert_eq!(String::from("TGCA*ACGT"), fasta_records[1].seq);
    }
}
