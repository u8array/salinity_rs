// std-Variante: vollständige Fehler mit thiserror und optional serde_json-Quellen
#[cfg(feature = "std")]
use thiserror::Error;

#[cfg(feature = "std")]
#[derive(Error, Debug)]
pub enum AppError {
    #[cfg(feature = "cli")]
    #[error("Error reading from stdin: {source}")]
    ReadStdin {
        #[source]
        source: std::io::Error,
    },

    #[cfg(feature = "cli")]
    #[error("Error reading file '{path}': {source}")]
    ReadFile {
        path: String,
        #[source]
        source: std::io::Error,
    },

    #[cfg(feature = "cli")]
    #[error("Invalid JSON for --inputs-json: {source}")]
    ParseInputsJson {
        #[source]
        source: serde_json::Error,
    },

    #[cfg(feature = "cli")]
    #[error("Invalid JSON for --assumptions-json: {source}")]
    ParseAssumptionsJson {
        #[source]
        source: serde_json::Error,
    },

    #[cfg(feature = "cli")]
    #[error("Invalid JSON in input document: {source}")]
    ParseCmdInputJson {
        #[source]
        source: serde_json::Error,
    },

    #[cfg(feature = "cli")]
    #[error("Could not serialize output to JSON: {source}")]
    SerializeOutput {
        #[source]
        source: serde_json::Error,
    },

    #[error("Unexpected error: {0}")]
    Other(String),

    #[cfg(feature = "cli")]
    #[error("Missing input data: provide --input or --inputs-json")]
    MissingInputData,

    #[cfg(feature = "cli")]
    #[error(
        "Missing assumptions: provide --assumptions-json or include 'assumptions' in the input document"
    )]
    MissingAssumptions,
}

#[cfg(not(feature = "std"))]
#[derive(Debug)]
pub enum AppError {
    Other(&'static str),
}
