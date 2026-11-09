export interface RealFunction {
  f: (x: number) => number;
  fprime: (x: number) => number;
}

/** f(x) = a(x − p)² + q, f′(x) = 2a(x − p). q defaults to 0. */
export function quad(a: number, p: number, q = 0): RealFunction {
  return {
    f: (x) => a * (x - p) ** 2 + q,
    fprime: (x) => 2 * a * (x - p),
  };
}

/** f(x) = x⁴ − 2x², f′(x) = 4x³ − 4x. Valleys at x = ±1. */
export function quartic(): RealFunction {
  return {
    f: (x) => x ** 4 - 2 * x ** 2,
    fprime: (x) => 4 * x ** 3 - 4 * x,
  };
}
