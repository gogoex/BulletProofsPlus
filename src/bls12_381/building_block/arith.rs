use std::sync::Once;
use mcl_rust::*;

pub struct Arith();

static INIT: Once = Once::new();

impl Arith {
  pub fn init() {
    INIT.call_once(|| {
      let b = init(CurveType::BLS12_381);
      if !b {
          panic!("Initializing mcl-rust library failed");
      }
    });
  }
}
