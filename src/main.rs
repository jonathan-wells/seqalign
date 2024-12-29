use std::io::Error;
use clap::Parser;

use seqalign::{Aligner, parse_fastafile};

#[derive(Parser)]
#[command(version, about)]
struct Args {
    #[arg(short, long)]
    query: String,

    #[arg(short, long)]
    target: String,

    #[arg(short, long)]
    scoring_matrix: String,

    #[arg(short='g', long)]
    open_penalty: i32,

    #[arg(short='x', long)]
    extension_penalty: i32,

    #[arg(short, long)]
    output: String,

    #[arg(short='m', long)]
    score: bool,

}

fn main() -> Result<(), Error> {
    let args = Args::parse();

    let aligner = Aligner::new(
        &args.scoring_matrix,
        args.open_penalty,
        args.extension_penalty
        )?;
    
    let queries = parse_fastafile(&args.query)?;
    let targets = parse_fastafile(&args.target)?;
    
    // if arecords.len() > 1 || brecords.len() >  1 {
    //     eprintln!("WARNING: seqalign currently only accepts single pairs.\n\
    //               Defaulting to first sequence per file.")
    // }
    
    for query_record in &queries {
        for target_record in &targets {
            let aresult = aligner.align(&query_record.seq, &target_record.seq);
            let alignment = aresult.alignment();
            println!(">{}\n{}", query_record.name, alignment.0);
            println!(">{}\n{}\n--\n", target_record.name, alignment.1);
            
            if args.score {
                println!("Score: {}", aresult.score());
            }
        }
    }
    
    Ok(())
}
