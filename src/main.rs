use std::io::Error;
use clap::Parser;

use seqalign::{Aligner, parse_fastafile};

#[derive(Parser, Debug)]
#[command(version, about)]
struct Args {
    #[arg(short, long)]
    afile: String,

    #[arg(short, long)]
    bfile: String,

    #[arg(short, long)]
    scoring_matrix: String,

    #[arg(short='o', long)]
    open_penalty: i32,

    #[arg(short, long)]
    extension_penalty: i32,

    #[arg(short='t', long)]
    output: String,

}

fn main() -> Result<(), Error> {
    let args = Args::parse();

    let aligner = Aligner::new(
        &args.scoring_matrix,
        args.open_penalty,
        args.extension_penalty
        )?;
    
    let arecords = parse_fastafile(&args.afile)?;
    let brecords = parse_fastafile(&args.bfile)?;
    
    if arecords.len() > 1 || brecords.len() >  1 {
        eprintln!("WARNING: seqalign currently only accepts single pairs.\n\
                  Defaulting to first sequence per file.")
    }
    // TODO: extend to make this accept multi-sequence files.
    let afile = &arecords[0];
    let bfile = &brecords[0];
    let aresult = aligner.align(&afile.seq, &bfile.seq);
    
    let alignment = aresult.alignment();
    println!(">{}\n{}", afile.name, alignment.0);
    println!(">{}\n{}", bfile.name, alignment.1);
    
    Ok(())
}
