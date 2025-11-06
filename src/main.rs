fn main() {
    if let Err(e) = salinity_rs::adapters::run() {
        eprintln!("{}", e);
        std::process::exit(1);
    }
}
