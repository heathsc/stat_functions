pub mod beta;
pub mod error;
pub mod pnorm;
pub mod students_t;

pub type Result<T> = std::result::Result<T, error::StatFuncError>;
