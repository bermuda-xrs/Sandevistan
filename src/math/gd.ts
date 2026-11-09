/** |x| above this, or a non-finite value, is divergence. Iteration stops. */
export const DIVERGENCE_LIMIT = 1e6;

export function isDivergent(value: number): boolean {
  return !Number.isFinite(value) || Math.abs(value) > DIVERGENCE_LIMIT;
}

/** One gradient step: x − η · f′(x). */
export function step(
  x: number,
  fprime: (x: number) => number,
  eta: number,
): number {
  return x - eta * fprime(x);
}

/**
 * [x0, x1, …] for up to n steps.
 * Stops after the first divergent value so later terms are not computed.
 */
export function run(
  x0: number,
  fprime: (x: number) => number,
  eta: number,
  n: number,
): number[] {
  const xs = [x0];
  if (n <= 0 || isDivergent(x0)) return xs;

  let x = x0;
  for (let i = 0; i < n; i++) {
    x = step(x, fprime, eta);
    xs.push(x);
    if (isDivergent(x)) break;
  }
  return xs;
}

/** Convergence ratio r = 1 − 2aη for f(x) = a(x − p)² + q. */
export function ratio(a: number, eta: number): number {
  return 1 - 2 * a * eta;
}

/** Closed form xₙ = p + rⁿ (x₀ − p). */
export function closedForm(x0: number, p: number, r: number, n: number): number {
  return p + r ** n * (x0 - p);
}
