use clap::{command, Parser, Subcommand, Args};

mod dataset;

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct CliParser {
    #[command(subcommand)]
    subcommand: CliSubcommand,
}

#[derive(Subcommand, Debug, Clone)]
enum CliSubcommand {
    /// Main entry point for the benchmark runner
    Bench(BenchArgs),

    /// Runs an individual algorithm on a data set
    RunJob, // TODO: runner args
}

#[derive(Args, Debug, Clone)]
struct BenchArgs {

}


fn bench(bench_args: &BenchArgs) {

}


fn main() {
    let args = CliParser::parse();

    match args.subcommand {
        CliSubcommand::Bench(ref bench_args) => bench(bench_args),
        CliSubcommand::RunJob => eprintln!("Not implemented yet!")
    }
}
