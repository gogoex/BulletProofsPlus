pub trait Zero<T> {
  fn zero(&self) -> T;
  fn is_zero(&self) -> bool;
}

pub trait Zero2<T> {
  fn zero() -> T;
  fn is_zero(&self) -> bool;
}
