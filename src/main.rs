use std::io::Error;
use clap::Parser;

use seqalign::{Aligner, parse_fastafile};

#[derive(Parser, Debug)]
#[command(version, about)]
struct Args {
    #[arg(short, long)]
    asequence: String,

    #[arg(short, long)]
    bsequence: String,

    #[arg(short, long)]
    scoring_matrix: String,

    #[arg(short, long)]
    open_penalty: i32,

    #[arg(short, long)]
    extension_penalty: i32,
}

fn main() -> Result<(), Error> {
    let args = Args::parse();

    let aligner = Aligner::new(
        args.scoring_matrix.as_str(),
        args.open_penalty,
        args.extension_penalty
        )?;
    
    let arecords = parse_fastafile(&args.asequence)?;
    let brecords = parse_fastafile(&args.bsequence)?;
    let asequence = &arecords[0];
    let bsequence = &brecords[0];
    let aresult = aligner.align(&asequence.seq, &bsequence.seq);
    
    let alignment = aresult.alignment();
    println!(">{}\n{}", asequence.name, alignment.0);
    println!(">{}\n{}", bsequence.name, alignment.1);
    
    Ok(())
}
