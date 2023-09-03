use thiserror::Error;

#[derive(Error, Debug)]
pub enum StatFuncError {
    #[error("Invalid parameters for beta function (alpha and beta must be >0")]
    InvalidBetaParameters,
    #[error("Invalid probability parameter (must be between 0 and 1")]
    InvalidProbability,
    #[error("Invalid degrees of freedom parameter for Students-t distribution (must be > 0)")]
    InvalidStudentsTParameter,
}
