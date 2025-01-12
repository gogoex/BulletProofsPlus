use crate::{
  impl_mul,
  impl_affine_add,
  secp256k1::building_block::{
    field::prime_field_elem::PrimeFieldElem,
    secp256k1::{
      jacobian_point::JacobianPoint,
      secp256k1::Secp256k1,
    },
    zero::Zero,
  },
};
use std::{
  fmt,
  ops::{Add, Mul},
  rc::Rc,
};

#[derive(Clone)]
pub struct AffinePoint {
  pub curve: Rc<Secp256k1>,
  pub x: PrimeFieldElem,
  pub y: PrimeFieldElem,
  pub is_inf: bool,
}

impl fmt::Debug for AffinePoint {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{{ x: {:?}, y: {:?}, is_inf: {:?} }}", &self.x, &self.y, self.is_inf)
  }
}

impl AffinePoint {
  pub fn new(curve: &Rc<Secp256k1>, x: &PrimeFieldElem, y: &PrimeFieldElem) -> Self {
    AffinePoint {
      curve: curve.clone(),
      x: x.clone(),
      y: y.clone(),
      is_inf: false,
    }
  }

  pub fn rand_point(&self, exclude_zero: bool) -> Self {
    let g = &self.curve.g();
    loop {
      let multiplier = self.curve.f_n.rand_elem(exclude_zero);
      let p = g * multiplier;
      if !exclude_zero || !p.is_zero() { return p; }
    }
  }

  pub fn inv(&self) -> Self {
    if self.is_inf {
      panic!("Cannot calculate the inverse of zero");
    }
    AffinePoint::new(
      &self.curve,
      &self.x,
      &self.y.inv(),
    )
  }
}

impl From<JacobianPoint> for AffinePoint {
  fn from(p: JacobianPoint) -> Self {
    if p.z.is_zero() {
      panic!("z is not expected to be zero");
    } else {
      let z2 = p.z.sq();
      let z3 = &z2 * &p.z;
      let x = &p.x / z2;
      let y = &p.y / z3;
      AffinePoint::new(
        &p.curve,
        &x,
        &y,
      )
    }
  }
}

impl Zero<AffinePoint> for AffinePoint {
  fn zero(&self) -> Self {
    AffinePoint {
      curve: self.curve.clone(),
      x: self.curve.f.elem(&0u8),
      y: self.curve.f.elem(&0u8),
      is_inf: true,
    }
  }

  fn is_zero(&self) -> bool {
    self.is_inf
  }
}

impl_affine_add!(AffinePoint);

impl_mul!(PrimeFieldElem, AffinePoint);

impl PartialEq for AffinePoint {
  fn eq(&self, other: &Self) -> bool {
    if self.is_inf != other.is_inf {  // false if one is zero and the other is non-zero
      false
    } else if self.is_inf {  // true if both are zero
      true
    } else {  // otherwise check if coordinates are the same
      self.x == other.x && self.y == other.y
    }
  }
}

impl Eq for AffinePoint {}

#[cfg(test)]
mod tests {
  use super::*;
  use num_bigint::BigUint;

  #[test]
  fn scalar_mul() {
    let curve = Rc::new(Secp256k1::new());
    let g = &curve.g();

    {
      let act = g * curve.f.elem(&1u8);
      assert_eq!(&act, g);
    }
    {
      let act = g * curve.f.elem(&2u8);
      let exp = g + g;
      assert_eq!(act, exp);
    }
    {
      let act = g * curve.f.elem(&3u8);
      let exp = g + g + g;
      assert_eq!(act, exp);
    }
  }

  #[test]
  fn add_same_point() {
    let curve = &Secp256k1::new();
    let g = &curve.g();
    let g2 = g + g;
    let exp_x = BigUint::parse_bytes(b"89565891926547004231252920425935692360644145829622209833684329913297188986597", 10).unwrap();
    let exp_y = BigUint::parse_bytes(b"12158399299693830322967808612713398636155367887041628176798871954788371653930", 10).unwrap();
    assert_eq!(g2.x.e, exp_x);
    assert_eq!(g2.y.e, exp_y);
  }

  #[test]
  fn add_same_point_y_eq_0() {
    // TODO implement this. need to find the x-coord when y is zero
  }

  #[test]
  fn add_vertical_line() {
    let curve = &Rc::new(Secp256k1::new());
    let g = curve.g();
    let a = g.clone();
    let b = AffinePoint::new(curve, &a.x, &-&a.y);
    let exp = g.zero();
    let act = &a + &b;
    assert_eq!(act, exp);
  }

  #[test]
  fn add_inf_and_affine() {
    let curve = &Secp256k1::new();
    let g = &curve.g();
    let inf = g.zero();
    let inf_plus_g = g + &inf;
    assert_eq!(g, &inf_plus_g);
  }

  #[test]
  fn add_affine_and_inf() {
    let curve = &Secp256k1::new();
    let g = &curve.g();
    let inf = g.zero();
    let g_plus_inf = &inf + g;
    assert_eq!(g, &g_plus_inf);
  }

  #[test]
  fn add_inf_and_inf() {
    let curve = Secp256k1::new();
    let g = &curve.g();
    let inf = g.zero();
    let inf_plus_inf = &inf + &inf;
    assert_eq!(inf_plus_inf, inf);
  }

  struct Xy<'a> {
    curve: &'a Rc<Secp256k1>,
    _n: &'a str,
    x: &'a [u8; 64],
    y: &'a [u8; 64],
  }

  impl<'a> Into<AffinePoint> for Xy<'a> {
    fn into(self) -> AffinePoint {
      let gx = BigUint::parse_bytes(self.x, 16).unwrap();
      let gy = BigUint::parse_bytes(self.y, 16).unwrap();
      AffinePoint::new(
        &self.curve,
        &PrimeFieldElem::new(&Rc::new(self.curve.f.clone()), &gx),
        &PrimeFieldElem::new(&Rc::new(self.curve.f.clone()), &gy),
      )
    }
  }

  // expects a + b = c
  struct AddTestCase {
    a: usize,
    b: usize,
    c: usize,
  }

  impl AddTestCase {
    fn new(a: usize, b: usize, c: usize) -> AddTestCase {
      if a + b != c {
        panic!("Bad add test case: {} + {} = {}", a, b, c);
      }
      AddTestCase { a, b, c }
    }
  }

  fn get_g_multiples<'a>(curve: &Rc<Secp256k1>) -> Vec<AffinePoint> {
    let ps = vec![
      Xy { curve, _n: "1", x: b"79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", y: b"483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8" },
      Xy { curve, _n: "2", x: b"C6047F9441ED7D6D3045406E95C07CD85C778E4B8CEF3CA7ABAC09B95C709EE5", y: b"1AE168FEA63DC339A3C58419466CEAEEF7F632653266D0E1236431A950CFE52A" },
      Xy { curve, _n: "3", x: b"F9308A019258C31049344F85F89D5229B531C845836F99B08601F113BCE036F9", y: b"388F7B0F632DE8140FE337E62A37F3566500A99934C2231B6CB9FD7584B8E672" },
      Xy { curve, _n: "4", x: b"E493DBF1C10D80F3581E4904930B1404CC6C13900EE0758474FA94ABE8C4CD13", y: b"51ED993EA0D455B75642E2098EA51448D967AE33BFBDFE40CFE97BDC47739922" },
      Xy { curve, _n: "5", x: b"2F8BDE4D1A07209355B4A7250A5C5128E88B84BDDC619AB7CBA8D569B240EFE4", y: b"D8AC222636E5E3D6D4DBA9DDA6C9C426F788271BAB0D6840DCA87D3AA6AC62D6" },
      Xy { curve, _n: "6", x: b"FFF97BD5755EEEA420453A14355235D382F6472F8568A18B2F057A1460297556", y: b"AE12777AACFBB620F3BE96017F45C560DE80F0F6518FE4A03C870C36B075F297" },
      Xy { curve, _n: "7", x: b"5CBDF0646E5DB4EAA398F365F2EA7A0E3D419B7E0330E39CE92BDDEDCAC4F9BC", y: b"6AEBCA40BA255960A3178D6D861A54DBA813D0B813FDE7B5A5082628087264DA" },
      Xy { curve, _n: "8", x: b"2F01E5E15CCA351DAFF3843FB70F3C2F0A1BDD05E5AF888A67784EF3E10A2A01", y: b"5C4DA8A741539949293D082A132D13B4C2E213D6BA5B7617B5DA2CB76CBDE904" },
      Xy { curve, _n: "9", x: b"ACD484E2F0C7F65309AD178A9F559ABDE09796974C57E714C35F110DFC27CCBE", y: b"CC338921B0A7D9FD64380971763B61E9ADD888A4375F8E0F05CC262AC64F9C37" },
      Xy { curve, _n: "10", x: b"A0434D9E47F3C86235477C7B1AE6AE5D3442D49B1943C2B752A68E2A47E247C7", y: b"893ABA425419BC27A3B6C7E693A24C696F794C2ED877A1593CBEE53B037368D7" },
    ];
    let g = &curve.g();
    let mut gs: Vec<AffinePoint> = vec![g.clone()];  // gs[0] is used to match index and g's n and will not be actually used
    for p in ps {
      gs.push(p.into());
    }
    gs
  }

  #[test]
  fn scalar_mul_smaller_nums() {
    let curve = Rc::new(Secp256k1::new());
    let g = &curve.g();
    let gs = get_g_multiples(&curve);

    for n in 1usize..=10 {
      let res = g * curve.f.elem(&n);
      assert_eq!(&res, &gs[n]);
    }
  }

  struct ScalarMulTest<'a> {
    k: &'a [u8; 64],
    x: &'a [u8; 64],
    y: &'a [u8; 64],
  }

  #[test]
  fn scalar_mul_gen_pubkey() {
    let test_cases = vec![
      ScalarMulTest {
        k: b"AA5E28D6A97A2479A65527F7290311A3624D4CC0FA1578598EE3C2613BF99522",
        x: b"34F9460F0E4F08393D192B3C5133A6BA099AA0AD9FD54EBCCFACDFA239FF49C6",
        y: b"0B71EA9BD730FD8923F6D25A7A91E7DD7728A960686CB5A901BB419E0F2CA232",
      },
      ScalarMulTest {
        k: b"7E2B897B8CEBC6361663AD410835639826D590F393D90A9538881735256DFAE3",
        x: b"D74BF844B0862475103D96A611CF2D898447E288D34B360BC885CB8CE7C00575",
        y: b"131C670D414C4546B88AC3FF664611B1C38CEB1C21D76369D7A7A0969D61D97D",
      },
      ScalarMulTest {
        k: b"6461E6DF0FE7DFD05329F41BF771B86578143D4DD1F7866FB4CA7E97C5FA945D",
        x: b"E8AECC370AEDD953483719A116711963CE201AC3EB21D3F3257BB48668C6A72F",
        y: b"C25CAF2F0EBA1DDB2F0F3F47866299EF907867B7D27E95B3873BF98397B24EE1",
      },
      ScalarMulTest {
        k: b"376A3A2CDCD12581EFFF13EE4AD44C4044B8A0524C42422A7E1E181E4DEECCEC",
        x: b"14890E61FCD4B0BD92E5B36C81372CA6FED471EF3AA60A3E415EE4FE987DABA1",
        y: b"297B858D9F752AB42D3BCA67EE0EB6DCD1C2B7B0DBE23397E66ADC272263F982",
      },
      ScalarMulTest {
        k: b"1B22644A7BE026548810C378D0B2994EEFA6D2B9881803CB02CEFF865287D1B9",
        x: b"F73C65EAD01C5126F28F442D087689BFA08E12763E0CEC1D35B01751FD735ED3",
        y: b"F449A8376906482A84ED01479BD18882B919C140D638307F0C0934BA12590BDE",
      },
    ];

    use std::time::Instant;
    let curve = Rc::new(Secp256k1::new());
    let g = &curve.g();

    for t in &test_cases {
      let k = BigUint::parse_bytes(t.k, 16).unwrap();
      let x = BigUint::parse_bytes(t.x, 16).unwrap();
      let y = BigUint::parse_bytes(t.y, 16).unwrap();
      let p = AffinePoint::new(
        &curve,
        &PrimeFieldElem::new(&Rc::new(curve.f.clone()), &x),
        &PrimeFieldElem::new(&Rc::new(curve.f.clone()), &y),
      );

      let beg = Instant::now();
      let gk = g * curve.f.elem(&k);
      let end = beg.elapsed();
      println!("Large number scalar mul done in {}.{:03} sec", end.as_secs(), end.subsec_nanos() / 1_000_000);
      assert_eq!(p, gk);
    }
  }

  #[test]
  fn add_different_points() {
    let curve = Rc::new(Secp256k1::new());
    let large_1 = Xy {
      curve: &curve,
      _n: "28948022309329048855892746252171976963209391069768726095651290785379540373584",
      x: b"A6B594B38FB3E77C6EDF78161FADE2041F4E09FD8497DB776E546C41567FEB3C",
      y: b"71444009192228730CD8237A490FEBA2AFE3D27D7CC1136BC97E439D13330D55",
    };
    let large_2 = Xy {
      curve: &curve,
      _n: "57896044618658097711785492504343953926418782139537452191302581570759080747168",
      x: b"00000000000000000000003B78CE563F89A0ED9414F5AA28AD0D96D6795F9C63",
      y: b"3F3979BF72AE8202983DC989AEC7F2FF2ED91BDD69CE02FC0700CA100E59DDF3",
    };
    let large_3 = Xy {
      curve: &curve,
      _n: "86844066927987146567678238756515930889628173209306178286953872356138621120752",
      x: b"E24CE4BEEE294AA6350FAA67512B99D388693AE4E7F53D19882A6EA169FC1CE1",
      y: b"8B71E83545FC2B5872589F99D948C03108D36797C4DE363EBD3FF6A9E1A95B10",
    };

    let gs = get_g_multiples(&curve);

    let test_cases = [
      AddTestCase::new(1, 2, 3),
      AddTestCase::new(2, 2, 4),
      AddTestCase::new(2, 6, 8),
      AddTestCase::new(3, 4, 7),
      AddTestCase::new(5, 1, 6),
      AddTestCase::new(5, 2, 7),
      AddTestCase::new(8, 1, 9),
      AddTestCase::new(9, 1, 10),
    ];

    for tc in test_cases {
      let res = &gs[tc.a] + &gs[tc.b];
      assert_eq!(&res, &gs[tc.c]);
    }

    let l1: AffinePoint = large_1.into();
    let l2: AffinePoint = large_2.into();
    let l3: AffinePoint = large_3.into();

    let l1_plus_l2 = &l1 + &l2;
    assert_eq!(l1_plus_l2, l3);
  }
}