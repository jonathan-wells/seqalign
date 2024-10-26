use std::io::Error;
use clap::Parser;

use seqalign::{Aligner, RognesAligner, parse_fastafile};

// #[derive(Parser)]
// #[command(version, about)]
// struct Args {
//     #[arg(short, long)]
//     query: String,

//     #[arg(short, long)]
//     target: String,

//     #[arg(short, long)]
//     scoring_matrix: String,

//     #[arg(short='g', long)]
//     open_penalty: i16,

//     #[arg(short='x', long)]
//     extension_penalty: i16,

//     #[arg(short, long)]
//     output: String,

// }

// fn main() -> Result<(), Error> {
//     let args = Args::parse();

//     let aligner = RognesAligner::new(
//         &args.scoring_matrix,
//         args.open_penalty,
//         args.extension_penalty
//         )?;
    
//     let queries = parse_fastafile(&args.query)?;
//     let targets = parse_fastafile(&args.target)?;
    
//     for query_record in &queries {
//         for target_record in &targets {
//             let aresult = aligner.align(&query_record.seq, &target_record.seq);
//             let alignment = aresult.alignment();
//             println!(">{}\n{}", query_record.name, alignment.0);
//             println!(">{}\n{}\n--\n", target_record.name, alignment.1);
//         }
//     }
    
//     Ok(())
// }

fn main() -> Result<(), Error> {

    let aligner = RognesAligner::new(
        "BLOSUM62",
        10,
        1
        )?;

    let query = "HEAG";
    let target = "HEAGPAW";
    let aresult = aligner.align(query, target);

    Ok(())
}
