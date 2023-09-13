use super::{
    beta::{betain, lbeta},
    Result,
};
use libm::lgamma;

use crate::error::StatFuncError;

pub struct StudentsT {
    v: f64,
    konst: f64,  // 1/(v*lnbeta)
    lnbeta: f64, // Beta(v/2, 1.2)
}

impl StudentsT {
    pub fn new(v: f64) -> Result<Self> {
        check_students_t_param(v)?;
        let lnbeta = lbeta(0.5 * v, 0.5).unwrap();
        let konst = -(0.5 * v.ln() + lnbeta);
        Ok(Self { v, konst, lnbeta })
    }

    pub fn dt(&self, t: f64) -> f64 {
        _students_t_pdf(t, self.v, self.konst)
    }

    pub fn ldt(&self, t: f64) -> f64 {
        _lstudents_t_pdf(t, self.v, self.konst)
    }

    pub fn pt(&self, t: f64) -> f64 {
        _students_t_cdf(t, self.v, Some(self.lnbeta))
    }
}

fn check_students_t_param(v: f64) -> Result<()> {
    if v <= 0.0 {
        Err(StatFuncError::InvalidStudentsTParameter)
    } else {
        Ok(())
    }
}

#[inline]
fn _lstudents_t_pdf(t: f64, v: f64, konst: f64) -> f64 {
    konst + (1.0 + t * t / v).ln() * -0.5 * (v + 1.0)
}

#[inline]
fn _students_t_pdf(t: f64, v: f64, konst: f64) -> f64 {
    _lstudents_t_pdf(t, v, konst).exp()
}

#[inline]
fn _students_t_cdf(t: f64, v: f64, lnbeta: Option<f64>) -> f64 {
    let x = v / (v + t * t);
    let z = 0.5 * betain(0.5 * v, 0.5, x, lnbeta).unwrap();
    if t < 0.0 {
        z
    } else {
        1.0 - z
    }
}

fn ldt(t: f64, v: f64) -> Result<f64> {
    check_students_t_param(v)?;
    let konst = lgamma(0.5 * (v + 1.0)) - lgamma(0.5 * v) - 0.5 * (v * std::f64::consts::PI).ln();
    Ok(_lstudents_t_pdf(t, v, konst))
}

fn dt(t: f64, v: f64) -> Result<f64> {
    ldt(t, v).map(|z| z.exp())
}

fn pt(t: f64, v: f64) -> Result<f64> {
    check_students_t_param(v)?;
    Ok(_students_t_cdf(t, v, None))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ldt_test() {
        let z = ldt(-3.0, 4.0).unwrap();
        assert!((z + 3.9274667438658417).abs() < 1.0e-12)
    }

    #[test]
    fn dt_test() {
        let z = dt(2.5, 2.8).unwrap();
        assert!((z - 0.039339657969571784).abs() < 1.0e-12)
    }

    #[test]
    fn ldt_test1() {
        let s = StudentsT::new(4.0).unwrap();
        let z = s.ldt(-3.0);
        assert!((z + 3.9274667438658417).abs() < 1.0e-12)
    }

    #[test]
    fn dt_test1() {
        let s = StudentsT::new(2.8).unwrap();
        let z = s.dt(2.5);
        assert!((z - 0.039339657969571784).abs() < 1.0e-12)
    }

    #[test]
    fn dt_test2() {
        let s = StudentsT::new(4.0).unwrap();
        let z = s.dt(0.0);
        assert!((z - 0.375).abs() < 1.0e-12)
    }

    #[test]
    fn pt_test() {
        let z = pt(2.5, 2.8).unwrap();
        assert!((z - 0.953134106244337).abs() < 1.0e-12)
    }

    #[test]
    fn pt_test1() {
        let s = StudentsT::new(2.8).unwrap();
        let z = s.pt(-3.4);
        println!("{z}");
        assert!((z - 0.02355410567174815).abs() < 1.0e-12)
    }
}
