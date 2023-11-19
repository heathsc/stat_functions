pub mod beta;
pub mod error;
mod pnorm;
pub mod students_t;

pub type Result<T> = std::result::Result<T, error::StatFuncError>;
