use phf::{phf_map};

// Define static strings for all implemented BLOSUM scoring matrices. 
pub static BLOSUM30: &str = include_str!("../resources/BLOSUM30");
pub static BLOSUM35: &str = include_str!("../resources/BLOSUM35");
pub static BLOSUM40: &str = include_str!("../resources/BLOSUM40");
pub static BLOSUM45: &str = include_str!("../resources/BLOSUM45");
pub static BLOSUM50: &str = include_str!("../resources/BLOSUM50");
pub static BLOSUM55: &str = include_str!("../resources/BLOSUM55");
pub static BLOSUM60: &str = include_str!("../resources/BLOSUM60");
pub static BLOSUM62: &str = include_str!("../resources/BLOSUM62");
pub static BLOSUM65: &str = include_str!("../resources/BLOSUM65");
pub static BLOSUM70: &str = include_str!("../resources/BLOSUM70");
pub static BLOSUM75: &str = include_str!("../resources/BLOSUM75");
pub static BLOSUM80: &str = include_str!("../resources/BLOSUM80");
pub static BLOSUM85: &str = include_str!("../resources/BLOSUM85");
pub static BLOSUM90: &str = include_str!("../resources/BLOSUM90");
pub static BLOSUM100: &str = include_str!("../resources/BLOSUM100");

pub const LANES: usize = 8;
pub const AA_ALPHABET_SIZE: usize = 24;
pub const AA_ALPHABET: [char; AA_ALPHABET_SIZE] = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'
];

pub static AA_ALPHABET_MAP: phf::Map<char, usize> = phf_map! {
    'A' => 0,
    'R' => 1,
    'N' => 2,
    'D' => 3,
    'C' => 4,
    'Q' => 5,
    'E' => 6,
    'G' => 7,
    'H' => 8,
    'I' => 9,
    'L' => 10,
    'K' => 11,
    'M' => 12,
    'F' => 13,
    'P' => 14,
    'S' => 15,
    'T' => 16,
    'W' => 17,
    'Y' => 18,
    'V' => 19,
    'B' => 20,
    'Z' => 21,
    'X' => 22,
    '*' => 23,
};
