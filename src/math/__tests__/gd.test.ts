import { describe, expect, it } from "vitest";
import {
  closedForm,
  isDivergent,
  quad,
  quartic,
  ratio,
  run,
} from "../index.ts";

function expectTrace(actual: number[], expected: number[]) {
  expect(actual).toHaveLength(expected.length);
  expected.forEach((value, index) => {
    expect(actual[index]).toBeCloseTo(value, 6);
  });
}

/** Deterministic [0, 1) generator so the closed-form sample stays reproducible. */
function mulberry32(seed: number): () => number {
  let state = seed >>> 0;
  return () => {
    state = (state + 0x6d2b79f5) >>> 0;
    let t = Math.imul(state ^ (state >>> 15), state | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), state | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

describe("4.4 검증값", () => {
  const valley = quad(1, 3);

  it("ratio(1, η) 는 0.5, −1.2, −1, 0", () => {
    expect(ratio(1, 0.25)).toBeCloseTo(0.5, 6);
    expect(ratio(1, 1.1)).toBeCloseTo(-1.2, 6);
    expect(ratio(1, 1)).toBeCloseTo(-1, 6);
    expect(ratio(1, 0.5)).toBeCloseTo(0, 6);
  });

  it("1. η = 0.25, x0 = 0 은 절반씩 3에 다가간다", () => {
    expectTrace(run(0, valley.fprime, 0.25, 6), [
      0, 1.5, 2.25, 2.625, 2.8125, 2.90625, 2.953125,
    ]);
  });

  it("2. η = 1.1, x0 = 0 은 커지며 부호가 바뀐다", () => {
    expectTrace(run(0, valley.fprime, 1.1, 4), [
      0, 6.6, -1.32, 8.184, -3.2208,
    ]);
  });

  it("3. η = 1/a 이면 r = −1 이라 0과 6을 왕복한다", () => {
    expectTrace(run(0, valley.fprime, 1, 4), [0, 6, 0, 6, 0]);
  });

  it("4. η = 1/(2a) 이면 한 걸음에 바닥에 도착한다", () => {
    expectTrace(run(0, valley.fprime, 0.5, 3), [0, 3, 3, 3]);
  });

  it("5. 임의의 이차함수에서 run 은 closedForm 과 모든 걸음이 일치한다", () => {
    const rng = mulberry32(0x5a17c0de);
    const samples = 120;
    const steps = 20;

    for (let sample = 0; sample < samples; sample++) {
      const a = (1 - rng()) * 5;
      const p = rng() * 40 - 20;
      const x0 = rng() * 40 - 20;
      const eta = rng() * 4;
      const r = ratio(a, eta);
      const xs = run(x0, quad(a, p).fprime, eta, steps);

      expect(a).toBeGreaterThan(0);
      expect(a).toBeLessThanOrEqual(5);

      for (let n = 0; n < xs.length; n++) {
        const actual = xs[n] ?? Number.NaN;
        const expected = closedForm(x0, p, r, n);
        const label = `sample ${sample} n=${n} a=${a} p=${p} x0=${x0} eta=${eta}`;
        if (Number.isFinite(actual) && Number.isFinite(expected)) {
          const scale = Math.max(1, Math.abs(expected));
          if (scale < 10) {
            expect(actual, label).toBeCloseTo(expected, 6);
          } else {
            expect(Math.abs(actual - expected) / scale, label).toBeLessThan(1e-8);
          }
        } else {
          expect(isDivergent(actual), label).toBe(true);
          expect(isDivergent(expected), label).toBe(true);
        }
      }

      const firstBad = xs.findIndex((value) => isDivergent(value));
      if (firstBad === -1) {
        expect(xs).toHaveLength(steps + 1);
      } else {
        expect(firstBad).toBe(xs.length - 1);
      }
    }
  });

  it("6. quartic 는 x0 = ±0.4, η = 0.1, 200걸음 뒤 ±1 에 가깝다", () => {
    const { fprime } = quartic();
    const right = run(0.4, fprime, 0.1, 200);
    const left = run(-0.4, fprime, 0.1, 200);
    const rightEnd = right[right.length - 1] ?? Number.NaN;
    const leftEnd = left[left.length - 1] ?? Number.NaN;

    expect(right).toHaveLength(201);
    expect(left).toHaveLength(201);
    expect(Math.abs(rightEnd - 1)).toBeLessThan(1e-6);
    expect(Math.abs(leftEnd - -1)).toBeLessThan(1e-6);
    expect(rightEnd).toBeCloseTo(1, 6);
    expect(leftEnd).toBeCloseTo(-1, 6);
  });

  it("7. quartic 는 x0 = 0 에서 모든 걸음이 0이다", () => {
    const xs = run(0, quartic().fprime, 0.1, 20);
    expect(xs).toHaveLength(21);
    for (const value of xs) expect(value).toBeCloseTo(0, 6);
  });

  it("8. quartic η = 0.2, x0 = 2 는 발산으로 판정된다", () => {
    const xs = run(2, quartic().fprime, 0.2, 30);
    const last = xs[xs.length - 1] ?? Number.NaN;

    expect(xs[0]).toBeCloseTo(2, 6);
    expect(xs[1]).toBeCloseTo(-2.8, 6);
    expect(xs.length).toBeLessThan(31);
    expect(isDivergent(last)).toBe(true);
    expect(xs.findIndex((value) => isDivergent(value))).toBe(xs.length - 1);
  });
});

describe("4.3 발산 처리", () => {
  it("절댓값이 1e6을 넘거나 NaN/Infinity 이면 발산이다", () => {
    expect(isDivergent(1e6)).toBe(false);
    expect(isDivergent(-1e6)).toBe(false);
    expect(isDivergent(1e6 + 1)).toBe(true);
    expect(isDivergent(-1e6 - 1)).toBe(true);
    expect(isDivergent(Number.NaN)).toBe(true);
    expect(isDivergent(Number.POSITIVE_INFINITY)).toBe(true);
    expect(isDivergent(Number.NEGATIVE_INFINITY)).toBe(true);
  });

  it("발산한 값 뒤로는 반복하지 않는다", () => {
    const xs = run(2, quartic().fprime, 0.2, 30);
    expect(xs.length).toBeGreaterThan(1);
    expect(isDivergent(xs[xs.length - 1] ?? Number.NaN)).toBe(true);
    expect(xs.slice(0, -1).every((value) => !isDivergent(value))).toBe(true);
  });
});
