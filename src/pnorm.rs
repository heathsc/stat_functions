use libm::ldexp;

pub fn pnorm(z: f64, lower_tail: bool) -> f64 {
    if z.is_infinite() {
        match (z.is_sign_negative(), lower_tail) {
            (true, true) => 0.0,   // -inf, lower tail
            (true, false) => 1.0,  // -inf, upper tail
            (false, true) => 1.0,  // +inf, lower tail
            (false, false) => 0.0, // +inf, upper tail
        }
    } else if z.is_nan() {
        z
    } else {
        _pnorm(z, lower_tail)
    }
}

const M_SQRT_32: f64 = 5.656_854_249_492_381;
const M_1_SQRT_2PI: f64 = 0.398_942_280_401_432_7;
fn _pnorm(x: f64, lower_tail: bool) -> f64 {
    const A: f64 = 0.065682337918207449113;
    const AB: [(f64, f64); 4] = [
        (2.2352520354606839287, 47.20258190468824187),
        (161.02823106855587881, 976.09855173777669322),
        (1067.6894854603709582, 10260.932208618978205),
        (18154.981253343561249, 45507.789335026729956),
    ];
    const C: f64 = 1.0765576773720192317e-8;
    const CD: [(f64, f64); 8] = [
        (0.39894151208813466764, 22.266688044328115691),
        (8.8831497943883759412, 235.38790178262499861),
        (93.506656132177855979, 1519.377599407554805),
        (597.27027639480026226, 6485.558298266760755),
        (2494.5375852903726711, 18615.571640885098091),
        (6848.1904505362823326, 34900.952721145977266),
        (11602.651437647350124, 38912.003286093271411),
        (9842.7148383839780218, 19685.429676859990727),
    ];
    const P: f64 = 0.02307344176494017303;
    const PQ: [(f64, f64); 5] = [
        (0.21589853405795699, 1.28426009614491121),
        (0.1274011611602473639, 0.468238212480865118),
        (0.022235277870649807, 0.0659881378689285515),
        (0.001421619193227893466, 0.00378239633202758244),
        (2.9112874951168792e-5, 7.29751555083966205e-5),
    ];
    const EPS: f64 = f64::EPSILON * 0.5;

    let do_del = |x, temp| -> f64 {
        let xsq = unsafe { ldexp(ldexp(x, 4).trunc(), -4) };
        let del = (x - xsq) * (x + xsq);
        (-xsq * xsq / 2.0 - del / 2.0).exp() * temp
    };

    let swap_tail = |x, p: f64| match (x < 0.0, lower_tail) {
        (true, true) => p,
        (true, false) => 1.0 - p,
        (false, true) => 1.0 - p,
        (false, false) => p,
    };
    let y = x.abs();
    if y <= 0.67448975 {
        let (xnum, xden) = if y > EPS {
            let xsq = x.powi(2);
            AB[..3].iter().fold((A * xsq, xsq), |(num, den), (a, b)| {
                ((num + a) * xsq, (den + b) * xsq)
            })
        } else {
            (0.0, 0.0)
        };
        let (a, b) = unsafe { AB.get_unchecked(3) };
        let temp = x * (xnum + a) / (xden + b);
        if lower_tail {
            0.5 + temp
        } else {
            0.5 - temp
        }
    } else if y <= M_SQRT_32 {
        let (xnum, xden) = CD[..7].iter().fold((C * y, y), |(num, den), (c, d)| {
            ((num + c) * y, (den + d) * y)
        });
        let (c, d) = unsafe { CD.get_unchecked(7) };
        let temp = (xnum + c) / (xden + d);
        swap_tail(x, do_del(y, temp))
    } else if (lower_tail && (-37.5193..8.2924).contains(&x))
        || (!lower_tail && (-8.2924..37.5193).contains(&x))
    {
        let xsq = x.powi(-2);
        let (xnum, xden) = PQ[..4].iter().fold((P * xsq, xsq), |(num, den), (p, q)| {
            ((num + p) * xsq, (den + q) * xsq)
        });
        let (p, q) = unsafe { PQ.get_unchecked(4) };
        let temp = (M_1_SQRT_2PI - xsq * (xnum + p) / (xden + q)) / y;
        swap_tail(x, do_del(x, temp))
    } else {
        swap_tail(x, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() {
        let p = pnorm(0.25, true);
        assert!((p - 0.5987063256829237).abs() < 1.0e-12)
    }
    #[test]
    fn test2() {
        let p = pnorm(-0.125, true);
        eprintln!("{p}");
        assert!((p - 0.45026177516988714).abs() < 1.0e-12)
    }
    #[test]
    fn test3() {
        let p = pnorm(-0.125, false);
        eprintln!("{p}");
        assert!((p - 0.5497382248301129).abs() < 1.0e-12)
    }
    #[test]
    fn test4() {
        let p = pnorm(1.96, true);
        eprintln!("{p}");
        assert!((p - 0.9750021048517796).abs() < 1.0e-12)
    }
    #[test]
    fn test5() {
        let p = pnorm(-3.0, false);
        eprintln!("{p}");
        assert!((p - 0.9986501019683699).abs() < 1.0e-12)
    }
    #[test]
    fn test6() {
        let p = pnorm(-25.0, true);
        eprintln!("{:e}", p.ln());
        assert!((p.ln() - -3.1663940800802027e2).abs() < 1.0e-12)
    }
    #[test]
    fn test7() {
        let p = pnorm(25.0, false);
        eprintln!("{:e}", p.ln());
        assert!((p.ln() - -3.1663940800802027e2).abs() < 1.0e-12)
    }

    #[test]
    fn test8() {
        let p = pnorm(25.0, true);
        eprintln!("{:e}", p);
        assert!(1.0 - p < 1.0e-12)
    }
}
