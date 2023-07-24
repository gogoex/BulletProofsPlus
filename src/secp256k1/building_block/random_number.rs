use rand_chacha::{ChaChaRng, rand_core::SeedableRng as SeedableRng2};

pub struct RandomNumber {
  pub gen: ChaChaRng,
}

impl RandomNumber {
  pub fn new() -> Self {
    let seed = [
        1, 0, 52, 0, 0, 0, 0, 0, 1, 0, 10, 0, 22, 32, 0, 0, 2, 0, 55, 49, 0, 11, 0, 0, 3, 0, 0, 0, 0,
        0, 2, 92,
    ];
    let gen = ChaChaRng::from_seed(seed);
    RandomNumber { gen }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use rand::RngCore;

  #[test]
  fn generate() {
    let mut r = RandomNumber::new();
    let mut buf = [0u8; 32];
    r.gen.fill_bytes(&mut buf);
  }
}