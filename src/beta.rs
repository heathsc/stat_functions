use super::{error::StatFuncError, Result};

use libm::lgamma;

#[inline]
fn _lbeta(p: f64, q: f64) -> f64 {
    lgamma(p) + lgamma(q) - lgamma(p + q)
}

pub fn lbeta(p: f64, q: f64) -> Result<f64> {
    check_beta_params(p, q)?;
    Ok(_lbeta(p, q))
}

pub fn beta(p: f64, q: f64) -> Result<f64> {
    lbeta(p, q).map(|x| x.exp())
}

/// Purpose:
///
///   betain() computes the incomplete Beta function ratio.
///
/// Licensing:
///
///   This code is distributed under the GNU LGPL license.
///
/// Modified:
///
///   31 October 2010
///
/// Author:
///
///   Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
///   C version by John Burkardt.
///   Rust version by Simon C Heath
/// Reference:
///
///   KL Majumder, GP Bhattacharjee,
///   Algorithm AS 63:
///   The incomplete Beta Integral,
///   Applied Statistics,
///   Volume 22, Number 3, 1973, pages 409-411.
///
/// Parameters:
///
///   Input, x, the argument, between 0 and 1.
///
///   Input, p, q, the parameters, which
///   must be positive.
///
///   Input, lnbeta, the logarithm of the complete beta function. If this
///   is None then the value is calculated.  Any supplied value is assumed
///   to be correct.
///
///   Output, the value of the incomplete Beta function ratio.
pub fn betain(p: f64, q: f64, x: f64, lnbeta: Option<f64>) -> Result<f64> {
    check_beta_params(p, q)?;
    if !(0.0..=1.0).contains(&x) {
        Err(StatFuncError::InvalidProbability)
    } else if x == 0.0 || x == 1.0 {
        Ok(x)
    } else {
        let lnbeta = lnbeta.unwrap_or_else(|| _lbeta(p, q));
        let accuracy = 1.0e-14;

        // Change tail if necessary
        let mut psq = p + q;
        let (xx, cx, pp, qq, flip) = if p < psq * x {
            (1.0 - x, x, q, p, true)
        } else {
            (x, 1.0 - x, p, q, false)
        };

        let (mut term, mut ai, mut value) = (1.0, 1.0, 1.0);
        let mut ns = (qq + cx * psq) as i64;

        // Use the Soper reduction formula
        let mut temp = qq - ai;
        let mut rx = if ns == 0 { xx } else { xx / cx };

        loop {
            term *= temp * rx / (pp + ai);
            value += term;
            temp = term.abs();

            if temp <= accuracy && temp <= accuracy * value {
                value *= (pp * xx.ln() + (qq - 1.0) * cx.ln() - lnbeta).exp() / pp;
                if flip {
                    value = 1.0 - value
                }
                break;
            }

            ai += 1.0;
            ns -= 1;
            match ns {
                1..=i64::MAX => temp = qq - ai,
                0 => {
                    temp = qq - ai;
                    rx = xx
                }
                _ => {
                    temp = psq;
                    psq += 1.0
                }
            }
        }
        Ok(value)
    }
}

fn check_beta_params(alpha: f64, beta: f64) -> Result<()> {
    if alpha <= 0.0 || beta <= 0.0 {
        Err(StatFuncError::InvalidBetaParameters)
    } else {
        Ok(())
    }
}

macro_rules! beta_inc {
    ($a:expr, $b:expr, $c:expr) => {
        betain($a, $b, $c, None)
    };
    ($a:expr, $b:expr, $c:expr, $beta:expr) => {
        betain($a, $b, $c, Some($beta))
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lnbeta_test() {
        let z = lbeta(4.0, 5.0).expect("Error in lbeta()");
        assert!((z + 5.63478960316925).abs() < 1.0e-12)
    }

    #[test]
    fn betain_test() {
        let z = betain(4.0, 5.0, 0.75, None).expect("Error in betain()");
        assert!((z - 0.9727020263671875).abs() < 1.0e-12)
    }

    #[test]
    fn betain_test1() {
        let z = betain(20.0, 5.0, 0.1, None).expect("Error in betain()");
        assert!((z - 7.1215255e-17).abs() < 1.0e-12)
    }

    #[test]
    fn betain_test2() {
        let z = betain(20.0, 5.0, 0.9, None).expect("Error in betain()");
        assert!((z - 0.914_925_114_121_329_2).abs() < 1.0e-12)
    }

    #[test]
    fn beta_inc() {
        let z = beta_inc!(20.0, 5.0, 0.9).expect("Error in betain()");
        assert!((z - 0.914_925_114_121_329_2).abs() < 1.0e-12)
    }

    #[test]
    fn beta_inc1() {
        let z = beta_inc!(20.0, 5.0, 0.9, lbeta(20.0, 5.0).unwrap()).expect("Error in betain()");
        assert!((z - 0.914_925_114_121_329_2).abs() < 1.0e-12)
    }
}
