# 재배열
def bit_reverse_permutation(a):
    n = len(a)
    bits = n.bit_length() - 1
    res = [0] * n

    for i in range(n):
        rev = 0
        x = i
        for _ in range(bits):
            rev = (rev << 1) | (x & 1)
            x >>= 1
        res[rev] = a[i]
    return res

# coeff -> value
def ntt_inplace(a, n, omega, q):
    a[:] = bit_reverse_permutation(a)
    length = 2
    while length <= n:
        half = length // 2
        w_m = pow(omega, n // length, q)
        for i in range(0, n, length):
            w = 1
            for j in range(half):
                u = a[i + j]
                v = (a[i + j + half] * w) % q
                a[i + j] = (u + v) % q
                a[i + j + half] = (u - v) % q
                w = (w * w_m) % q
        length *= 2

def intt_inplace(a, n, omega_inv, q):
    a[:] = bit_reverse_permutation(a)
    length = 2
    while length <= n:
        half = length // 2
        w_m = pow(omega_inv, n // length, q)
        for i in range(0, n, length):
            w = 1
            for j in range(half):
                u = a[i + j]
                v = a[i + j + half]
                t = (w * v) % q
                a[i + j] = (u + t) % q
                a[i + j + half] = (u - t) % q
                w = (w * w_m) % q
        length *= 2
    n_inv = pow(n, q - 2, q)
    for i in range(n):
        a[i] = (a[i] * n_inv) % q


class NTTPolynomialRing:
    def __init__(self, coefficients, n, q, psi):
        self.n = n
        self.q = q
        coeffs = [0] * n
        for i, c in enumerate(coefficients):
            coeffs[i % n] = (coeffs[i % n] + c) % q
        self.coefficients = [int(x % q) for x in coeffs]
        self.psi = psi
        self.psi_inv = pow(psi, q - 2, q)
        self.omega = pow(psi, 2, q)        # primitive n-th root
        self.omega_inv = pow(self.omega, q - 2, q)

    def copy(self):
        return NTTPolynomialRing(self.coefficients[:], self.n, self.q, self.psi)

    def _twist(self, arr):
        out = [0] * self.n
        cur = 1
        for i in range(self.n):
            out[i] = (arr[i] * cur) % self.q
            cur = (cur * self.psi) % self.q
        return out

    def _untwist(self, arr):
        out = [0] * self.n
        cur = 1
        psi_inv = self.psi_inv
        for i in range(self.n):
            out[i] = (arr[i] * cur) % self.q
            cur = (cur * psi_inv) % self.q
        return out

    def ntt(self):
        a = self._twist(self.coefficients)
        ntt_inplace(a, self.n, self.omega, self.q)
        return a

    def intt(self, values):
        b = values[:]
        intt_inplace(b, self.n, self.omega_inv, self.q)
        return self._untwist(b)

    def __mul__(self, other):
        A = self.ntt()
        B = other.ntt()
        C = [(A[i] * B[i]) % self.q for i in range(self.n)]
        c_coeffs = self.intt(C)
        return NTTPolynomialRing(c_coeffs, self.n, self.q, self.psi)

    def __add__(self, other):
        return NTTPolynomialRing([(self.coefficients[i] + other.coefficients[i]) % self.q for i in range(self.n)],
                                 self.n, self.q, self.psi)

    def coeffs(self):
        return self.coefficients[:]



# =========================
# NTT timing harness
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

    # ===== 파라미터 =====
    n = 2048
    q = 12289
    psi = 7   # 필요하면 바꿔

    # 랜덤 입력
    a_coeffs = [random.randrange(-q, q) for _ in range(n)]
    b_coeffs = [random.randrange(-q, q) for _ in range(n)]

    # NTT Ring
    a = NTTPolynomialRing(a_coeffs, n, q, psi)
    b = NTTPolynomialRing(b_coeffs, n, q, psi)

    print(f"Params: n={n}, q={q}")
    print(f"Nonzero count a={sum(1 for x in a.coefficients if x != 0)}, "
          f"b={sum(1 for x in b.coefficients if x != 0)}")

    # 1회 시간
    _, _ = time_once("NTT MUL (once)", lambda: a * b)

    # 반복 평균 시간
    print("\n[Repeat timing for stable average]")
    _, _, _ = time_repeat("NTT MUL (repeat)", lambda: a * b, repeats=50, warmup=10)