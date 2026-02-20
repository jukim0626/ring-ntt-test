class PolynomialRing:
    def __init__(self, coefficients, n, q):

        self.n = n
        self.q = q

        if len(coefficients) < n:
            self.coefficients = coefficients + [0] * (n - len(coefficients))
        elif len(coefficients) > n:
            self.coefficients = self._reduce_modulo(coefficients)
        else:
            self.coefficients = coefficients[:]

        # 계수 mod q
        self.coefficients = [c % self.q for c in self.coefficients]

    def _reduce_modulo(self, coeffs):

        result = [0] * self.n

        for i, coef in enumerate(coeffs):
            if i < self.n:
                result[i] += coef
            else:
                quotient = i // self.n
                remainder = i % self.n

                if quotient % 2 == 1:
                    result[remainder] -= coef
                else:
                    result[remainder] += coef

        return result

    # O(n)
    def __add__(self, other):
        if self.n != other.n or self.q != other.q:
            raise ValueError("Ring parameter error")

        result_coeffs = []
        for i in range(self.n):
            result_coeffs.append((self.coefficients[i] + other.coefficients[i]) % self.q)

        return PolynomialRing(result_coeffs, self.n, self.q)

    # O(n)
    def __sub__(self, other):
        if self.n != other.n or self.q != other.q:
            raise ValueError("Ring parameter error")

        result_coeffs = []
        for i in range(self.n):
            result_coeffs.append((self.coefficients[i] - other.coefficients[i]) % self.q)

        return PolynomialRing(result_coeffs, self.n, self.q)

    # O(n^2)
    def __mul__(self, other):
        if self.n != other.n or self.q != other.q:
            raise ValueError("Ring parameter error")

        temp_result = [0] * (2 * self.n - 1)

        for i in range(self.n):
            if self.coefficients[i] == 0:
                continue
            for j in range(self.n):
                if other.coefficients[j] == 0:
                    continue
                temp_result[i + j] += self.coefficients[i] * other.coefficients[j]


        result_coeffs = self._reduce_modulo(temp_result)

        return PolynomialRing(result_coeffs, self.n, self.q)



# =========================
# Quick timing harness
# =========================
import time
import random


def time_once(label, fn):
    t0 = time.perf_counter()
    out = fn()
    t1 = time.perf_counter()
    dt = t1 - t0
    print(f"{label}: {dt:.6f} s")
    return out, dt


def time_repeat(label, fn, repeats=50, warmup=5):
    # warmup
    for _ in range(warmup):
        fn()

    t0 = time.perf_counter()
    out = None
    for _ in range(repeats):
        out = fn()
    t1 = time.perf_counter()

    total = t1 - t0
    avg = total / repeats
    print(f"{label}: total {total:.6f} s / {repeats} runs  -> avg {avg*1e6:.2f} µs")
    return out, total, avg


if __name__ == "__main__":
    # ====== 네가 여기만 바꾸면 됨 ======
    n = 2048
    q = 12289

    # 임의 입력 (원하는 숫자로 바꿔서 테스트)

    a_coeffs = [random.randrange(-q, q) for _ in range(n)]
    b_coeffs = [random.randrange(-q, q) for _ in range(n)]

    # ==================================
    a = PolynomialRing(a_coeffs, n, q)
    b = PolynomialRing(b_coeffs, n, q)

    print(f"Params: n={n}, q={q}")
    print(f"Nonzero count a={sum(1 for x in a.coefficients if x != 0)}, "
          f"b={sum(1 for x in b.coefficients if x != 0)}")

    # 1회 시간
    _, _ = time_once("ADD (once)", lambda: a + b)
    _, _ = time_once("MUL (once, naive O(n^2))", lambda: a * b)

    # 반복 평균 시간(더 안정적)
    print("\n[Repeat timing for stable average]")
    _, _, _ = time_repeat("ADD (repeat)", lambda: a + b, repeats=200, warmup=20)
    _, _, _ = time_repeat("MUL (repeat, naive O(n^2))", lambda: a * b, repeats=30, warmup=5)


