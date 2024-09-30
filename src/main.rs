use std::io::Error;

use seqalign::Aligner;

fn main() -> Result<(), Error> {
    let aligner = Aligner::new("BLOSUM62", -5, 1)?;

    let aresult = aligner.align(String::from("HEAGHPAW"), String::from("HEWAGHQPAW"));

    let alignment = aresult.alignment();
    println!("Alignment");
    println!("{}", alignment.0);
    println!("{}\n", alignment.1);

    println!("Scores");
    for i in 0..aresult.s1.len() {
        for j in 0..aresult.s2.len() {
            print!("{:3} ", aresult.matrix[i][j]);
        }
        println!("");
    }

    println!("\nTraceback");
    for i in 0..aresult.s1.len() {
        for j in 0..aresult.s2.len() {
            print!("{} ", aresult.traceback_matrix[i][j]);
        }
        println!("");
    }
    Ok(())
}
