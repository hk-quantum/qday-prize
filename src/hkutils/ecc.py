# Copyright 2025 Hirokazu Kobayashi
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import Tuple, List
from qiskit.circuit import Qubit
from . import stdgates as g
from . import quantum as qu
from .quantum import QubitTarget, IGate, CompositeGate
import numpy as np
import sys
import math


def get_pow2_route(target: int) -> List[Tuple[int, int]]:
    route = []
    bit_num = target.bit_length()
    for i in range(bit_num - 1):
        route.append((i, i))
    for i in range(bit_num - 1):
        if target & (1 << i):
            route.append((i, len(route)))
    return route


def minus_x_to_y_mod_n(x: QubitTarget, y: QubitTarget, n: int) -> qu.IGate:
    """
    y: |0> -> |-x mod n>
    """
    qc = CompositeGate()
    bit_num = n.bit_length()
    # y = n
    qc <<= a_to_x(n, y)
    # y -= x
    qc <<= qu.inv @ x_add_z_to_z(x, y)
    # 0だったら n をリセット
    qc <<= qu.neg_ctrl(*x) @ a_to_x(n, y)
    return qc


def x_pow_2_mod_n_to_zc(x: QubitTarget, n: int, zc: QubitTarget) -> IGate:
    """
    z: |0> -> |x^2 mod n>
    """
    qc = CompositeGate()
    bit_num = n.bit_length()
    qc = qu.CompositeGate()
    if len(zc) < bit_num + 1:
        raise ValueError("zのビット数が足りません")
    # n未満であれば x_i^2 はそのまま設定できる
    cnt = 0
    total = 0
    for i in range(bit_num):
        total |= 1 << (i * 2)
        if total >= n:
            break
        cnt += 1
    for i in range(cnt):
        # そのまま設定
        qc <<= qu.ctrl(x[i]) @ g.x(zc[i * 2])
    for i in range(bit_num):
        for j in range(i, bit_num):
            if i == j:
                if i < cnt:
                    continue
                qc <<= qu.ctrl(x[i]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], ((1 << i) * (1 << j)) % n, n
                )
            else:
                qc <<= qu.ctrl(x[i], x[j]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], (2 * (1 << i) * (1 << j)) % n, n
                )
    return qc


def x_pow_2_add_zc_to_zc(x: QubitTarget, zc: QubitTarget, n: int) -> IGate:
    """
    |x>|y>|z> -> |x>|y>|(z + x^2) mod n>
    """
    bit_num = n.bit_length()
    qc = qu.CompositeGate()
    if len(zc) < bit_num + 1:
        raise ValueError("zのビット数が足りません")
    for i in range(bit_num):
        for j in range(i, bit_num):
            if i == j:
                qc <<= qu.ctrl(x[i]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], ((1 << i) * (1 << j)) % n, n
                )
            else:
                qc <<= qu.ctrl(x[i], x[j]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], (2 * (1 << i) * (1 << j)) % n, n
                )
    return qc


def x_pow_2_mul_a_mod_n_to_zc(x: QubitTarget, a: int, zc: QubitTarget, n: int) -> IGate:
    """
    |x>|y>|z> -> |x>|y>|(z + a*x^2) mod n>
    """
    qc = CompositeGate()
    bit_num = n.bit_length()
    qc = qu.CompositeGate()
    if len(zc) < bit_num + 1:
        raise ValueError("zのビット数が足りません")
    qc <<= qu.ctrl(x[0]) @ a_to_x(a, zc[:bit_num])
    for i in range(bit_num):
        for j in range(i, bit_num):
            if i == j:
                if i == 0:
                    continue
                qc <<= qu.ctrl(x[i]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], (a * (1 << i) * (1 << j)) % n, n
                )
            else:
                qc <<= qu.ctrl(x[i], x[j]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], (2 * a * (1 << i) * (1 << j)) % n, n
                )
    return qc


def x_pow_2_mul_a_add_zc_mod_n_to_zc(
    x: QubitTarget, a: int, zc: QubitTarget, n: int
) -> IGate:
    """
    |x>|y>|z> -> |x>|y>|(z + a*x^2) mod n>
    """
    qc = CompositeGate()
    bit_num = n.bit_length()
    qc = qu.CompositeGate()
    if len(zc) < bit_num + 1:
        raise ValueError("zのビット数が足りません")
    for i in range(bit_num):
        for j in range(i, bit_num):
            if i == j:
                qc <<= qu.ctrl(x[i]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], (a * (1 << i) * (1 << j)) % n, n
                )
            else:
                qc <<= qu.ctrl(x[i], x[j]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], (2 * a * (1 << i) * (1 << j)) % n, n
                )
    return qc


def x_pow_3_mod_n_to_zc(x: QubitTarget, n: int, zc: QubitTarget) -> IGate:
    """
    z: |0> -> |x^3 mod n>

    (a + b + c)^3 = a^3 + b^3 + c^3 + 3a^2b + 3a^2c + 3b^2a + 3b^2c + 3c^2a + 3c^2b + 6abc

    """
    qc = CompositeGate()
    bit_num = n.bit_length()
    qc = qu.CompositeGate()
    if len(zc) < bit_num + 1:
        raise ValueError("zのビット数が足りません")
    # 1の位はそのまま設定
    qc <<= qu.ctrl(x[0]) @ g.x(zc[0])
    for i in range(bit_num):
        for j in range(i, bit_num):
            mul = 3
            ctrl_regs = [x[i]]
            if j != i:
                # 3*a*b^2
                ctrl_regs.append(x[j])
                qc <<= qu.ctrl(*ctrl_regs) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], (3 * (1 << i) * (1 << j) * (1 << j)) % n, n
                )
                mul = 6
            elif i > 0:
                # a^3
                qc <<= qu.ctrl(x[i]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], ((1 << i) * (1 << j) * (1 << j)) % n, n
                )
            for k in range(j + 1, bit_num):
                qc <<= qu.ctrl(*ctrl_regs, x[k]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], (mul * (1 << i) * (1 << j) * (1 << k)) % n, n
                )
    return qc


def x_mul_a_mod_n_to_zc(x: QubitTarget, a: int, zc: QubitTarget, n: int) -> IGate:
    """
    z: |0> -> |(a*x) mod n>
    """
    bit_num = n.bit_length()
    if len(zc) < bit_num + 1:
        raise ValueError("z is too small " + str(len(zc)))
    qc = CompositeGate()
    d = a
    for i in range(bit_num):
        if i == 0:
            qc <<= qu.ctrl(x[i]) @ a_to_x(d, zc)
        else:
            qc <<= qu.ctrl(x[i]) @ xc_add_a_mod_n_to_xc(zc, d, n)
        d = (d << 1) % n
    return qc


def x_mul_a_add_zc_mod_n_to_zc(
    x: QubitTarget, a: int, zc: QubitTarget, n: int
) -> IGate:
    """
    |x>|y>|z> -> |x>|y>|(z + a*x) mod n>
    """
    bit_num = n.bit_length()
    if len(zc) < bit_num + 1:
        raise ValueError("z is too small " + str(len(zc)))
    qc = CompositeGate()
    d = a
    for i in range(bit_num):
        qc <<= qu.ctrl(x[i]) @ xc_add_a_mod_n_to_xc(zc, d, n)
        d = (d << 1) % n
    return qc


def x_mul_y_mod_n_to_zc(
    x: QubitTarget, y: QubitTarget, zc: QubitTarget, n: int
) -> qu.IGate:
    """
    z: |0> -> |(x*y) mod n>
    """
    bit_num = n.bit_length()
    if len(zc) < bit_num + 1:
        raise ValueError("zのビット数が足りません")
    qc = qu.CompositeGate()
    # 1の位はそのまま設定
    qc <<= qu.ctrl(x[0]) @ qu.ctrl(y[:bit_num]) @ g.x(zc[:bit_num])
    for i in range(1, bit_num):
        for j in range(bit_num):
            qc <<= qu.ctrl(x[i], y[j]) @ xc_add_a_mod_n_to_xc(
                zc[: bit_num + 1], ((1 << i) * (1 << j)) % n, n
            )
    return qc


def x_mul_y_add_zc_mod_n_to_zc(
    x: QubitTarget, y: QubitTarget, zc: QubitTarget, n: int
) -> qu.IGate:
    """
    |x>|y>|z> -> |x>|y>|(z + x*y) mod n>
    """
    bit_num = n.bit_length()
    if len(zc) < bit_num + 1:
        raise ValueError("zのビット数が足りません")
    qc = qu.CompositeGate()
    for i in range(bit_num):
        for j in range(bit_num):
            qc <<= qu.ctrl(x[i], y[j]) @ xc_add_a_mod_n_to_xc(
                zc[: bit_num + 1], ((1 << i) * (1 << j)) % n, n
            )
    return qc


def x_pow_a_mod_n_to_z_with_divs(
    x: QubitTarget,
    a: int,
    n: int,
    z: QubitTarget,
    anc: QubitTarget,
    divs: List[Tuple[int, int]],
) -> IGate:
    bit_num = n.bit_length()
    qc = CompositeGate()
    regs = [x]
    for i in range(1, len(divs)):
        regs.append(anc[i * bit_num : (i + 1) * bit_num])
    regs.append(z)
    carry = anc[len(divs) * bit_num : len(divs) * bit_num + 1]
    for i in range(len(divs)):
        dt = divs[i]
        if dt[0] == dt[1]:
            qc <<= x_pow_2_mod_n_to_zc(regs[dt[0]], n, regs[i + 1] + carry)
        else:
            qc <<= x_mul_y_mod_n_to_zc(regs[dt[0]], regs[dt[1]], regs[i + 1] + carry, n)
    return qc


def qft_dagger(reg: QubitTarget) -> qu.IGate:
    qc = CompositeGate()
    n = len(reg)
    for ix in range(n // 2):
        qc <<= g.swap(reg[ix], reg[n - ix - 1])
    for j in range(n):
        for m in range(j):
            qc <<= g.cp(-np.pi / float(2 ** (j - m)))(reg[m], reg[j])
        qc <<= g.h(reg[j])
    return qc


def x_inc_to_x(reg: QubitTarget) -> qu.IGate:
    qc = qu.CompositeGate()
    for i in reversed(range(1, len(reg))):
        qc <<= qu.ctrl(*reg[:i]) @ g.x(reg[i])
    qc <<= g.x(reg[0])
    return qc


def a_to_x(a: int, x: QubitTarget) -> IGate:
    qc = CompositeGate()
    bit_num = a.bit_length()
    for i in range(bit_num):
        if (a >> i) & 1:
            qc <<= g.x(x[i])
    return qc


def x_add_a_to_x(x: QubitTarget, a: int) -> IGate:
    qc = CompositeGate()
    bit_num = a.bit_length()
    # 量子ゲート数の見込みを計算
    num = len(x)
    pl_num = 0
    mi_num = 0
    b = (1 << len(x)) - a
    for i in range(len(x)):
        if a & (1 << i):
            pl_num += num
        if b & (1 << i):
            mi_num += num
        num -= 1
    if pl_num > mi_num:
        # 引き算の方がビット数が少ない
        for i in range(len(x)):
            if (b >> i) & 1:
                qc <<= qu.inv @ x_inc_to_x(x[i:])
    else:
        for i in reversed(range(bit_num)):
            if (a >> i) & 1:
                qc <<= x_inc_to_x(x[i:])
    return qc


def x_add_z_to_z(x: QubitTarget, z: QubitTarget) -> IGate:
    qc = CompositeGate()
    for i in range(min(len(x), len(z))):
        qc <<= qu.ctrl(x[i]) @ x_inc_to_x(z[i:])
    return qc


def xc_add_a_mod_n_to_xc(xc: QubitTarget, a: int, n: int) -> IGate:
    """
    xはbit_num+1
    """
    bit_num = n.bit_length()
    if len(xc) < bit_num + 1:
        raise ValueError("x is too small")
    qc = CompositeGate()
    # 1bitだけなら効率的にできる
    plus_bit = -1
    minus_bit = -1
    for i in range(bit_num):
        if (1 << i) == a:
            plus_bit = i
            break
        elif (1 << i) == n - a:
            minus_bit = i
            break
    if plus_bit >= 0:
        qc <<= x_add_a_to_x(xc[: bit_num + 1], a)
        qc <<= qu.inv @ x_add_a_to_x(xc[: bit_num + 1], n)
        qc <<= qu.ctrl(xc[bit_num]) @ x_add_a_to_x(xc[:bit_num], n)
        qc <<= qu.neg_ctrl(*xc[plus_bit:bit_num]) @ g.x(xc[bit_num])
        qc <<= g.x(xc[bit_num])
    elif minus_bit >= 0:
        qc <<= g.x(xc[bit_num])
        qc <<= qu.neg_ctrl(*xc[minus_bit:bit_num]) @ g.x(xc[bit_num])
        qc <<= qu.ctrl(xc[bit_num]) @ qu.inv @ x_add_a_to_x(xc[:bit_num], n)
        qc <<= x_add_a_to_x(xc[: bit_num + 1], n)
        qc <<= qu.inv @ x_add_a_to_x(xc[: bit_num + 1], n - a)
    else:
        qc <<= x_add_a_to_x(xc[: bit_num + 1], a)
        qc <<= qu.inv @ x_add_a_to_x(xc[: bit_num + 1], n)
        qc <<= qu.ctrl(xc[bit_num]) @ x_add_a_to_x(xc[:bit_num], n)
        qc <<= qu.inv @ x_add_a_to_x(xc[: bit_num + 1], a)
        qc <<= x_add_a_to_x(xc[:bit_num], a)
        qc <<= g.x(xc[bit_num])
    return qc


def minus_xc_mod_n_to_xc(xc: QubitTarget, n: int) -> IGate:
    qc = CompositeGate()
    bit_num = n.bit_length()
    if len(xc) < bit_num + 1:
        raise ValueError("x is too small")
    qc <<= g.x(xc[: bit_num + 1])
    qc <<= x_inc_to_x(xc[: bit_num + 1])
    qc <<= qu.ctrl(xc[bit_num]) @ x_add_a_to_x(xc[:bit_num], n)
    qc <<= qu.neg_ctrl(*xc[:bit_num]) @ g.x(xc[bit_num])
    qc <<= g.x(xc[bit_num])
    return qc


def x_add_zc_mod_n_to_zc(x: QubitTarget, zc: QubitTarget, n: int) -> IGate:
    bit_num = n.bit_length()
    if len(zc) < bit_num + 1:
        raise ValueError("z is too small " + str(len(zc)))
    qc = CompositeGate()
    qc <<= x_add_z_to_z(x[:bit_num], zc[: bit_num + 1])
    qc <<= qu.inv @ x_add_a_to_x(zc[: bit_num + 1], n)
    qc <<= qu.ctrl(zc[bit_num]) @ x_add_a_to_x(zc[:bit_num], n)
    qc <<= qu.inv @ x_add_z_to_z(x[:bit_num], zc[: bit_num + 1])
    qc <<= x_add_z_to_z(x[:bit_num], zc[:bit_num])
    qc <<= g.x(zc[bit_num])
    return qc


def cx_mul_2_mod_n_to_xc(xc: QubitTarget, n: int) -> IGate:
    qc = CompositeGate()
    bit_num = n.bit_length()
    if len(xc) < bit_num + 1:
        raise ValueError("x is too small")

    return qc


def x_pow_a_mod_n_to_z0(x: QubitTarget, a: int, n: int, z: QubitTarget) -> IGate:
    """
    zにゴミが残る
    x^(a & ~1)までの計算結果を z の下位ビットに入れる
    """
    bit_num = n.bit_length()
    a_bit_num = a.bit_length()
    qc = CompositeGate()
    if a == 0:
        # 0乗
        qc <<= g.x(z[0])
        return qc
    elif a == 1:
        qc <<= qu.ctrl(x[0:bit_num]) @ g.x(z[0:bit_num])
        return qc
    elif a == 2:
        qc <<= x_pow_2_mod_n_to_zc(x, n, z)
        return qc
    elif a == 3:
        qc <<= x_pow_3_mod_n_to_zc(x, n, z)
        return qc
    mul_list = [x[:bit_num]]
    # ２乗を繰り返す
    for i in range(1, a_bit_num):
        mul_list.append(z[(i - 1) * bit_num : bit_num * i])
        qc <<= x_pow_2_mod_n_to_zc(mul_list[-2], n, mul_list[-1] + z[bit_num * i :])
    tmpres = z[bit_num * (a_bit_num - 1) :]
    for i in reversed(range(1, a_bit_num - 1)):
        if a & (1 << i):
            # 掛け算する
            qc <<= x_mul_y_mod_n_to_zc(mul_list[i], mul_list[i + 1], tmpres, n)
            # 0に戻す
            qc <<= qu.inv @ x_pow_2_mod_n_to_zc(
                mul_list[i - 1], n, mul_list[i] + tmpres[bit_num:]
            )
            qc <<= g.swap(mul_list[i], tmpres[:bit_num])
        else:
            qc <<= g.swap(mul_list[i], mul_list[i + 1])
    if a & 1:
        qc <<= x_mul_y_mod_n_to_zc(x, mul_list[1], tmpres, n)
        qc <<= g.swap(mul_list[1], tmpres[:bit_num])
    return qc


def ctrl_with_state(reg: QubitTarget, state: int) -> qu.GateModifier:
    mod = None
    for i in range(len(reg)):
        if state & (1 << i):
            if mod is None:
                mod = qu.ctrl(reg[i])
            else:
                mod = mod @ qu.ctrl(reg[i])
        else:
            if mod is None:
                mod = qu.neg_ctrl(reg[i])
            else:
                mod = mod @ qu.neg_ctrl(reg[i])
    return mod  # type: ignore


class EccQuantum:
    a: int
    b: int
    G: Tuple[int, int]
    p: int
    order: int
    bit_num: int
    inv_size: int
    inv_list: List[Tuple[int, int]]

    def __init__(self, a: int, b: int, G: Tuple[int, int], p: int):
        self.a = a
        self.b = b
        self.G = G
        self.p = p
        self.bit_num = p.bit_length()
        self.order = 1
        P = self.G
        while P != (0, 0):
            P = self.add(P, self.G)
            self.order += 1
        self.inv_list = get_pow2_route(self.p - 2)
        self.inv_size = len(self.inv_list) + 1

    def get_nG(self, n: int) -> Tuple[int, int]:
        if n == 0:
            return (0, 0)
        elif n == 1:
            return self.G
        elif n < 0:
            x, y = self.get_nG(-n)
            return (x, -y % self.p)
        elif n % 2 == 0:
            x1, y1 = self.get_nG(n // 2)
            m = (3 * x1 * x1 + self.a) * pow(2 * y1, self.p - 2, self.p) % self.p
            x3 = (m * m - 2 * x1) % self.p
            y3 = (m * (x1 - x3) - y1) % self.p
            return (x3, y3)
        else:
            x1, y1 = self.get_nG(n - 1)
            x2, y2 = self.G
            if x1 == 0 and y1 == 0:
                return (x2, y2)
            if x2 == 0 and y2 == 0:
                return (x1, y1)
            if x1 == x2 and y1 == (-y2 % self.p):
                return (0, 0)
            m = (y2 - y1) * pow(x2 - x1, self.p - 2, self.p) % self.p
            x3 = (m * m - x1 - x2) % self.p
            y3 = (m * (x1 - x3) - y1) % self.p
            return (x3, y3)

    def get_aP(self, a: int, P: Tuple[int, int]) -> Tuple[int, int]:
        if a == 0:
            return (0, 0)
        elif a == 1:
            return P
        elif a < 0:
            x, y = self.get_aP(-a, P)
            return (x, -y % self.p)
        elif a % 2 == 0:
            x1, y1 = self.get_aP(a // 2, P)
            m = (3 * x1 * x1 + self.a) * pow(2 * y1, self.p - 2, self.p) % self.p
            x3 = (m * m - 2 * x1) % self.p
            y3 = (m * (x1 - x3) - y1) % self.p
            return (x3, y3)
        else:
            x1, y1 = self.get_aP(a - 1, P)
            x2, y2 = P
            if x1 == 0 and y1 == 0:
                return (x2, y2)
            if x2 == 0 and y2 == 0:
                return (x1, y1)
            if x1 == x2 and y1 == (-y2 % self.p):
                return (0, 0)
            m = (y2 - y1) * pow(x2 - x1, self.p - 2, self.p) % self.p
            x3 = (m * m - x1 - x2) % self.p
            y3 = (m * (x1 - x3) - y1) % self.p
            return (x3, y3)

    def add(self, P: Tuple[int, int], Q: Tuple[int, int]) -> Tuple[int, int]:
        if P == (0, 0):
            return Q
        if Q == (0, 0):
            return P
        if P[0] == Q[0] and P[1] == (-Q[1] % self.p):
            return (0, 0)
        if P == Q:
            m = (3 * P[0] * P[0] + self.a) * pow(2 * P[1], self.p - 2, self.p) % self.p
        else:
            m = (Q[1] - P[1]) * pow(Q[0] - P[0], self.p - 2, self.p) % self.p
        x3 = (m * m - P[0] - Q[0]) % self.p
        y3 = (m * (P[0] - x3) - P[1]) % self.p
        return (x3, y3)

    def is_valid(self, P: Tuple[int, int]) -> bool:
        if P == (0, 0):
            return True
        x, y = P
        return (y * y - (x * x * x + self.a * x + self.b)) % self.p == 0

    def __str__(self) -> str:
        return f"a={self.a} b={self.b} p={self.p} G={self.G} order={self.order}"

    def get_add_ancilla_size(self) -> int:
        """
        1回あたりの座標足し算で必要なAncillaビット
        これ以外に xx, yy, zz の bit_num * 6 が必要
        """
        return self.bit_num * (1 + self.inv_size) + 5

    def calc_xx_add_Q_lambda_dx_inv_to_z0(
        self, xx: QubitTarget, Q: Tuple[int, int], z: QubitTarget, carry: QubitTarget
    ) -> IGate:
        """
        lambdaの分母の逆元
        """
        qc = CompositeGate()
        dx = z[self.bit_num : self.bit_num * 2]
        qc <<= qu.ctrl(xx[: self.bit_num]) @ g.x(dx)
        qc <<= qu.inv @ xc_add_a_mod_n_to_xc(dx + carry, Q[0], self.p)
        qc <<= x_pow_a_mod_n_to_z0(
            dx[: self.bit_num],
            self.p - 2,
            self.p,
            z[: self.bit_num] + z[self.bit_num * 2 :] + carry,
        )
        return qc

    def calc_xx_add_Q_lambda_dy_to_yc(
        self, xx: QubitTarget, Q: Tuple[int, int], yc: QubitTarget
    ) -> IGate:
        """
        lambdaの分子 y2-y1
        """
        qc = CompositeGate()
        qc <<= qu.ctrl(xx[self.bit_num : self.bit_num * 2]) @ g.x(yc[: self.bit_num])
        qc <<= qu.inv @ xc_add_a_mod_n_to_xc(yc[: self.bit_num + 1], Q[1], self.p)
        return qc

    def calc_xx_add_Q_flag_to_y(
        self, xx: QubitTarget, Q: Tuple[int, int], y: QubitTarget
    ) -> IGate:
        """
        y[0] = x=Oフラグ
        y[1] = x+Q=Oフラグ
        """
        qc = CompositeGate()
        # x==O の無限遠点フラグ
        qc <<= qu.neg_ctrl(*xx[: self.bit_num * 2]) @ g.x(y[0])
        # x+Q==O の無限遠点フラグ
        qc <<= ctrl_with_state(
            xx[: self.bit_num * 2], ((self.p - Q[1]) << self.bit_num) | Q[0]
        ) @ g.x(y[1])
        return qc

    def calc_xx_add_yy_to_flag_lambda(
        self,
        xx: QubitTarget,
        yy: QubitTarget,
        flag: QubitTarget,
        lam: QubitTarget,
        ext: QubitTarget,
    ) -> IGate:
        """
        flag[0] = xx=Oフラグ
        flag[1] = yy=Oフラグ
        flag[2] = xx+yy=Oフラグ
        lam:
        ext:
          dx:
          dy:
          same: xx==yy フラグ
          carry: 1bit
          dx^-1:
          ancilla
        """
        qc = CompositeGate()
        # xx==O の無限遠点フラグ
        qc <<= qu.neg_ctrl(*xx[: self.bit_num * 2]) @ g.x(flag[0])
        # xy==O の無限遠点フラグ
        qc <<= qu.neg_ctrl(*yy[: self.bit_num * 2]) @ g.x(flag[1])

        # dx,dyの計算
        qc <<= qu.ctrl(xx[: self.bit_num * 2]) @ g.x(ext[: self.bit_num * 2])
        # dx
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(
            yy[: self.bit_num], ext[: self.bit_num] + ext[self.bit_num * 2 :], self.p
        )
        # dy
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(
            yy[self.bit_num : self.bit_num * 2],
            ext[self.bit_num :],
            self.p,
        )

        # x+Q==O の無限遠点フラグ(途中計算は戻す)
        # xx.y + yy.y = p の判定
        qc <<= qu.ctrl(xx[self.bit_num : self.bit_num * 2]) @ g.x(
            ext[self.bit_num * 2 : self.bit_num * 3]
        )
        qc <<= x_add_z_to_z(
            yy[self.bit_num : self.bit_num * 2],
            ext[self.bit_num * 2 : self.bit_num * 3],
        )
        qc <<= qu.inv @ x_add_a_to_x(ext[self.bit_num * 2 : self.bit_num * 3], self.p)
        # 判定
        qc <<= qu.neg_ctrl(
            *ext[: self.bit_num], *ext[self.bit_num * 2 : self.bit_num * 3]
        ) @ g.x(flag[2])
        # 戻し
        qc <<= x_add_a_to_x(ext[self.bit_num * 2 : self.bit_num * 3], self.p)
        qc <<= qu.inv @ x_add_z_to_z(
            yy[self.bit_num : self.bit_num * 2],
            ext[self.bit_num * 2 : self.bit_num * 3],
        )
        qc <<= (
            qu.inv
            @ qu.ctrl(xx[self.bit_num : self.bit_num * 2])
            @ g.x(ext[self.bit_num * 2 : self.bit_num * 3])
        )

        # 同一点の計算 dx=2*y1,dy=3x^2+a
        qc <<= qu.neg_ctrl(*ext[: self.bit_num * 2]) @ g.x(ext[self.bit_num * 2])
        # dx = 2*y1
        qc <<= qu.ctrl(ext[self.bit_num * 2]) @ x_mul_a_mod_n_to_zc(
            xx[self.bit_num : self.bit_num * 2],
            2,
            ext[: self.bit_num] + ext[self.bit_num * 2 + 1 :],
            self.p,
        )
        # dy=3*x^2+a
        qc <<= x_pow_2_mod_n_to_zc(
            xx[: self.bit_num],
            self.p,
            ext[self.bit_num * 2 + 1 :],
        )
        qc <<= qu.ctrl(ext[self.bit_num * 2]) @ x_mul_a_mod_n_to_zc(
            ext[self.bit_num * 2 + 1 : self.bit_num * 3 + 1],
            3,
            ext[self.bit_num : self.bit_num * 2] + ext[self.bit_num * 3 + 1 :],
            self.p,
        )
        if self.a > 0:
            qc <<= qu.ctrl(ext[self.bit_num * 2]) @ xc_add_a_mod_n_to_xc(
                ext[self.bit_num : self.bit_num * 2] + ext[self.bit_num * 3 + 1 :],
                self.a,
                self.p,
            )
        # 戻し
        qc <<= qu.inv @ x_pow_2_mod_n_to_zc(
            xx[: self.bit_num],
            self.p,
            ext[self.bit_num * 2 + 1 :],
        )

        # dx^-1の計算
        qc <<= x_pow_a_mod_n_to_z0(
            ext[: self.bit_num], self.p - 2, self.p, ext[self.bit_num * 2 + 2 :]
        )

        # lambdaの計算
        qc <<= x_mul_y_mod_n_to_zc(
            ext[self.bit_num * 2 + 2 : self.bit_num * 3 + 2],
            ext[self.bit_num : self.bit_num * 2],
            lam + ext[self.bit_num * 2 + 1 : self.bit_num * 2 + 2],
            self.p,
        )

        return qc

    def calc_xx_add_Q_lambda_to_z0(
        self, xx: QubitTarget, Q: Tuple[int, int], z: QubitTarget, carry: QubitTarget
    ) -> IGate:
        qc = CompositeGate()
        if len(z) < self.bit_num * 3 + 1:
            raise ValueError("z is too small")
        # 同じ点フラグ
        qc <<= ctrl_with_state(
            xx[: self.bit_num * 2], (Q[1] << self.bit_num) | Q[0]
        ) @ g.x(z[self.bit_num])
        # 同じ点の場合
        same_lambda = (
            (3 * Q[0] * Q[0] + self.a) * pow(2 * Q[1], self.p - 2, self.p) % self.p
        )
        qc <<= qu.ctrl(z[self.bit_num]) @ a_to_x(
            same_lambda,
            z[: self.bit_num],
        )
        # 違う点の場合
        qc <<= self.calc_xx_add_Q_lambda_dx_inv_to_z0(
            xx, Q, z[self.bit_num * 2 + 1 :], carry
        )
        qc <<= self.calc_xx_add_Q_lambda_dy_to_yc(
            xx, Q, z[self.bit_num + 1 : self.bit_num * 2 + 1] + carry
        )
        qc <<= qu.neg_ctrl(z[self.bit_num]) @ x_mul_y_mod_n_to_zc(
            z[self.bit_num * 2 + 1 : self.bit_num * 3 + 1],
            z[self.bit_num + 1 : self.bit_num * 2 + 1],
            z[: self.bit_num] + carry,
            self.p,
        )
        return qc

    def xx_add_Q_with_lambda_to_zz(
        self,
        xx: QubitTarget,
        Q: Tuple[int, int],
        lam: QubitTarget,
        zz: QubitTarget,
        ancilla: QubitTarget,
    ) -> IGate:
        qc = CompositeGate()

        # x3の計算
        # lambda^2
        qc <<= x_pow_2_mod_n_to_zc(lam, self.p, zz + ancilla)
        # -x1
        qc <<= qu.inv @ xc_add_a_mod_n_to_xc(zz[: self.bit_num + 1], Q[0], self.p)
        # -x2
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(xx[: self.bit_num], zz, self.p)

        # y3の計算
        # lambda * x1
        qc <<= x_mul_a_mod_n_to_zc(lam, Q[0], zz[self.bit_num :] + ancilla, self.p)
        # lambda * x3(戻し対象)
        qc <<= x_mul_y_mod_n_to_zc(
            lam,
            zz[: self.bit_num],
            ancilla,
            self.p,
        )

        # -lambda * x3
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(
            ancilla[: self.bit_num],
            zz[self.bit_num : self.bit_num * 2] + ancilla[self.bit_num :],
            self.p,
        )
        # -y1
        qc <<= qu.inv @ xc_add_a_mod_n_to_xc(
            zz[self.bit_num : self.bit_num * 2] + ancilla[self.bit_num :],
            Q[1],
            self.p,
        )

        # ここから戻し
        qc <<= qu.inv @ x_mul_y_mod_n_to_zc(
            lam,
            zz[: self.bit_num],
            ancilla,
            self.p,
        )

        return qc

    def xx_add_yy_with_lambda_to_zz(
        self,
        xx: QubitTarget,
        yy: QubitTarget,
        lam: QubitTarget,
        zz: QubitTarget,
        ancilla: QubitTarget,
    ) -> IGate:
        qc = CompositeGate()

        # x3の計算
        # lambda^2
        qc <<= x_pow_2_mod_n_to_zc(lam, self.p, zz + ancilla)
        # -x1
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(yy[: self.bit_num], zz, self.p)
        # -x2
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(xx[: self.bit_num], zz, self.p)

        # y3の計算
        # lambda * x1
        qc <<= x_mul_y_mod_n_to_zc(
            lam, xx[: self.bit_num], zz[self.bit_num :] + ancilla, self.p
        )
        # lambda * x3(戻し対象)
        qc <<= x_mul_y_mod_n_to_zc(
            lam,
            zz[: self.bit_num],
            ancilla,
            self.p,
        )

        # -lambda * x3
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(
            ancilla[: self.bit_num],
            zz[self.bit_num : self.bit_num * 2] + ancilla[self.bit_num :],
            self.p,
        )
        # -y1
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(
            xx[self.bit_num : self.bit_num * 2],
            zz[self.bit_num : self.bit_num * 2] + ancilla[self.bit_num :],
            self.p,
        )

        # ここから戻し
        qc <<= qu.inv @ x_mul_y_mod_n_to_zc(
            lam,
            zz[: self.bit_num],
            ancilla,
            self.p,
        )

        return qc

    def xx_add_Q_to_zz(
        self, xx: QubitTarget, Q: Tuple[int, int], zz: QubitTarget, ancilla: QubitTarget
    ) -> IGate:
        if len(xx) < self.bit_num * 2:
            raise ValueError("x is too small")
        qc = CompositeGate()
        carry = ancilla[0:1]
        flag = ancilla[1:3]
        lam = ancilla[self.bit_num + 3 : self.bit_num * 2 + 3]

        qc <<= self.calc_xx_add_Q_lambda_to_z0(
            xx, Q, ancilla[self.bit_num + 3 :], carry
        )
        qc <<= self.calc_xx_add_Q_flag_to_y(xx, Q, flag)
        # 最終座標
        qc <<= qu.neg_ctrl(*flag) @ self.xx_add_Q_with_lambda_to_zz(
            xx, Q, lam, zz[: self.bit_num * 2], ancilla[3 : self.bit_num + 3] + carry
        )
        # x=Oの時
        qc <<= qu.ctrl(flag[0]) @ a_to_x((Q[1] << self.bit_num) | Q[0], zz)

        # 戻し
        qc <<= qu.inv @ self.calc_xx_add_Q_flag_to_y(xx, Q, flag)
        qc <<= qu.inv @ self.calc_xx_add_Q_lambda_to_z0(
            xx, Q, ancilla[self.bit_num + 3 :], carry
        )
        return qc

    def xx_add_Q_to_zz0(
        self, xx: QubitTarget, Q: Tuple[int, int], zz: QubitTarget, ancilla: QubitTarget
    ) -> IGate:
        if len(xx) < self.bit_num * 2:
            raise ValueError("x is too small")
        # (x2, y2), lam, a1, carry
        # lam=(y2-y1)/(x2-x1)
        # x3=lam^2-x1-x2
        # y3=lam*(x1-x3)-y1
        # 手順
        # f1=(x1==0 && y1==0) (x1,y1)=O
        # x3=-x2
        # x2-=x1 (x2が使えなくなる)
        # y2-=y1 (y2が使えなくなる)
        # f2=(x2==0) (結果 O)
        # f3=(y2==0) 同一点
        # ctrl(f3) @ x(f2)
        # a1=x2^-1
        # lam=y2*a1
        # lam= (3*x1^2+a)/(2*y1) (if x1==x2)
        # x3+=lam^2
        # x3-=x1
        # y3=-y1
        # y3+=lam*x1
        # y3-=lam*x3
        # x3^=-x1^x1(f1=1)
        # y3^=-y1^y1(f1=1)
        qc = CompositeGate()
        x2, y2 = xx[: self.bit_num], xx[self.bit_num : self.bit_num * 2]
        x1, y1 = Q[0], Q[1]
        x3, y3 = zz[: self.bit_num], zz[self.bit_num : self.bit_num * 2]
        carry = ancilla[0:1]
        f1 = ancilla[1:2]
        f2 = ancilla[2:3]
        f3 = ancilla[3:4]
        lam = ancilla[4 : self.bit_num + 4]
        a1 = ancilla[self.bit_num + 4 :]

        # f1=(x1==0 && y1==0) (x1,y1)=O
        qc <<= qu.neg_ctrl(*x2, *y2) @ g.x(f1)
        # x3=-x2
        qc <<= minus_x_to_y_mod_n(x2, x3, self.p)
        # x2-=x1 (x2が使えなくなる)
        qc <<= xc_add_a_mod_n_to_xc(x2 + carry, self.p - x1, self.p)
        # y2-=y1 (y2が使えなくなる)
        qc <<= xc_add_a_mod_n_to_xc(y2 + carry, self.p - y1, self.p)
        # f2=(x2==0) (結果 O)
        qc <<= qu.neg_ctrl(*x2) @ g.x(f2)
        # f3=(y2==0) 同一点
        qc <<= qu.neg_ctrl(*y2) @ g.x(f3)
        # ctrl(f3) @ x(f2)
        qc <<= qu.ctrl(f3) @ g.x(f2)
        # a1=x2^-1
        qc <<= x_pow_a_mod_n_to_z_with_divs(
            x2,
            self.p - 2,
            self.p,
            a1[: self.bit_num],
            a1[self.bit_num :] + carry,
            self.inv_list,
        )
        # lam=y2*a1
        qc <<= qu.neg_ctrl(f1) @ x_mul_y_mod_n_to_zc(
            y2, a1[: self.bit_num], lam + carry, self.p
        )
        # lam= (3*x1^2+a)/(2*y1) (if x1==x2)
        lam_same = (3 * x1 * x1 + self.a) * pow(2 * y1, self.p - 2, self.p) % self.p
        qc <<= qu.ctrl(f3) @ a_to_x(lam_same, lam)
        # x3+=lam^2
        qc <<= x_pow_2_add_zc_to_zc(lam, x3 + carry, self.p)
        # x3-=x1
        qc <<= qu.inv @ xc_add_a_mod_n_to_xc(x3 + carry, x1, self.p)
        # y3=-y1
        qc <<= a_to_x(self.p - y1, y3)
        # y3+=lam*x1
        qc <<= x_mul_a_add_zc_mod_n_to_zc(lam, x1, y3 + carry, self.p)
        # y3-=lam*x3
        qc <<= qu.inv @ x_mul_y_add_zc_mod_n_to_zc(lam, x3, y3 + carry, self.p)
        # x3^=-x1^x1(f1=1)
        qc <<= qu.ctrl(f1) @ a_to_x((self.p - x1) ^ x1, x3)
        # y3^=-y1^y1(f1=1)
        qc <<= qu.ctrl(f1) @ a_to_x((self.p - y1) ^ y1, y3)
        # x3=0(f2=1)
        qc <<= qu.ctrl(f2) @ g.swap(x2, x3)
        # y3^=-x1(f2=1)
        qc <<= qu.ctrl(f2) @ a_to_x(self.p - y1, y3)

        return qc

    def _xx_add_yy_to_zz0(
        self, xx: QubitTarget, yy: QubitTarget, zz: QubitTarget, ancilla: QubitTarget
    ) -> IGate:
        qc = CompositeGate()
        x1, y1 = xx[: self.bit_num], xx[self.bit_num : self.bit_num * 2]
        x2, y2 = yy[: self.bit_num], yy[self.bit_num : self.bit_num * 2]
        x3, y3 = zz[: self.bit_num], zz[self.bit_num : self.bit_num * 2]
        f1 = ancilla[0:1]
        f2 = ancilla[1:2]
        f3 = ancilla[2:3]
        f4 = ancilla[3:4]
        carry = ancilla[4:5]
        lam = ancilla[5 : self.bit_num + 5]
        aa = ancilla[self.bit_num + 5 :]

        # f1=(x1==0 && y1==0) (x1,y1)=O
        qc <<= qu.neg_ctrl(*x1, *y1) @ g.x(f1)
        # f2=(x2==0 && y2==0) (x2,y2)=O
        qc <<= qu.neg_ctrl(*x2, *y2) @ g.x(f2)
        # ctrl(f1) x3=x2
        qc <<= qu.ctrl(f1, x2) @ g.x(x3)
        # ctrl(f1) y3=y2
        qc <<= qu.ctrl(f1, y2) @ g.x(y3)
        # ctrl(f2) x3=x1
        qc <<= qu.ctrl(f2, x1) @ g.x(x3)
        # ctrl(f2) y3=y1
        qc <<= qu.ctrl(f2, y1) @ g.x(y3)
        # negctrl(f1) x3=-x2
        qc <<= qu.neg_ctrl(f1) @ minus_x_to_y_mod_n(x2, x3, self.p)
        # x2-=x1 (mod p)
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(x1, x2 + carry, self.p)
        # y2-=y1 (mod p)
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(
            y1,
            y2 + carry,
            self.p,
        )
        # f3=(x2==0 && y2!=0) 結果 O
        # f4=(x2==0 && y2==0) 2倍
        qc <<= qu.neg_ctrl(*x2) @ g.x(f3)
        qc <<= qu.ctrl(f3) @ qu.neg_ctrl(*y2) @ g.x(f4)
        qc <<= qu.ctrl(f4) @ g.x(f3)
        # ctrl(f3) x3 = -x1 (inv)
        qc <<= qu.ctrl(f3) @ qu.inv @ minus_x_to_y_mod_n(x1, x3, self.p)
        # ctrl(f4) y2= (3*x1^2) + a (mod p)
        if self.a > 0:
            qc <<= qu.ctrl(f4) @ a_to_x(self.a, y2)
            qc <<= qu.ctrl(f4) @ x_pow_2_mul_a_add_zc_mod_n_to_zc(
                x1, 3, y2 + carry, self.p
            )
        else:
            qc <<= qu.ctrl(f4) @ x_pow_2_mul_a_mod_n_to_zc(x1, 3, y2 + carry, self.p)
        # ctrl(f4) x2= 2*y1 (mod p)
        qc <<= qu.ctrl(f4) @ x_mul_a_mod_n_to_zc(y1, 2, x2 + carry, self.p)
        # a2=x2^-1 (mod p)
        qc <<= x_pow_a_mod_n_to_z_with_divs(
            x2,
            self.p - 2,
            self.p,
            aa[: self.bit_num],
            aa[self.bit_num :] + carry,
            self.inv_list,
        )
        # lambda=y2*a2 (mod p)
        qc <<= x_mul_y_mod_n_to_zc(y2, aa, lam + carry, self.p)
        # negctrl(f1, f2, f3) x3 += lambda^2 - x1 (mod p)
        qc <<= qu.neg_ctrl(f1, f2, f3) @ x_pow_2_add_zc_to_zc(lam, x3 + carry, self.p)
        qc <<= (
            qu.neg_ctrl(f1, f2, f3)
            @ qu.inv
            @ x_add_zc_mod_n_to_zc(x1, x3 + carry, self.p)
        )
        # x1 -= x3 (mod p)
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(x3, x1 + carry, self.p)
        # negctrl(f1, f2, f3) y3 = lambda * x1 - y1 (mod p)
        qc <<= qu.neg_ctrl(f1, f2, f3) @ minus_x_to_y_mod_n(y1, y3, self.p)
        qc <<= qu.neg_ctrl(f1, f2, f3) @ x_mul_y_add_zc_mod_n_to_zc(
            lam, x1, y3 + carry, self.p
        )
        return qc

    def _xx_add_yy_to_zz0_no_dbl(
        self, xx: QubitTarget, yy: QubitTarget, zz: QubitTarget, ancilla: QubitTarget
    ) -> IGate:
        """
        doubling無視バージョン
        """
        qc = CompositeGate()
        x1, y1 = xx[: self.bit_num], xx[self.bit_num : self.bit_num * 2]
        x2, y2 = yy[: self.bit_num], yy[self.bit_num : self.bit_num * 2]
        x3, y3 = zz[: self.bit_num], zz[self.bit_num : self.bit_num * 2]
        f1 = ancilla[0:1]
        f2 = ancilla[1:2]
        f3 = ancilla[2:3]
        carry = ancilla[3:4]
        lam = ancilla[4 : self.bit_num + 4]
        aa = ancilla[self.bit_num + 4 :]

        # f1=(x1==0 && y1==0) (x1,y1)=O
        qc <<= qu.neg_ctrl(*x1, *y1) @ g.x(f1)
        # f2=(x2==0 && y2==0) (x2,y2)=O
        qc <<= qu.neg_ctrl(*x2, *y2) @ g.x(f2)
        # ctrl(f1) x3=x2
        qc <<= qu.ctrl(f1, x2) @ g.x(x3)
        # ctrl(f1) y3=y2
        qc <<= qu.ctrl(f1, y2) @ g.x(y3)
        # ctrl(f2) x3=x1
        qc <<= qu.ctrl(f2, x1) @ g.x(x3)
        # ctrl(f2) y3=y1
        qc <<= qu.ctrl(f2, y1) @ g.x(y3)
        # negctrl(f1) x3=-x2
        qc <<= qu.neg_ctrl(f1) @ minus_x_to_y_mod_n(x2, x3, self.p)
        # x2-=x1 (mod p)
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(x1, x2 + carry, self.p)
        # y2-=y1 (mod p)
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(
            y1,
            y2 + carry,
            self.p,
        )
        # f3=(x2==0) 結果 O doublingなし
        qc <<= qu.neg_ctrl(*x2) @ g.x(f3)
        # ctrl(f3) x3 = -x1 (inv)
        qc <<= qu.ctrl(f3) @ qu.inv @ minus_x_to_y_mod_n(x1, x3, self.p)
        # aa=x2^-1 (mod p)
        qc <<= x_pow_a_mod_n_to_z_with_divs(
            x2,
            self.p - 2,
            self.p,
            aa[: self.bit_num],
            aa[self.bit_num :] + carry,
            self.inv_list,
        )
        # lambda=y2*aa (mod p)
        qc <<= x_mul_y_mod_n_to_zc(y2, aa, lam + carry, self.p)
        # negctrl(f1, f2, f3) x3 += lambda^2 - x1 (mod p)
        qc <<= qu.neg_ctrl(f1, f2, f3) @ x_pow_2_add_zc_to_zc(lam, x3 + carry, self.p)
        qc <<= (
            qu.neg_ctrl(f1, f2, f3)
            @ qu.inv
            @ x_add_zc_mod_n_to_zc(x1, x3 + carry, self.p)
        )
        # x1 -= x3 (mod p)
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(x3, x1 + carry, self.p)
        # negctrl(f1, f2, f3) y3 = lambda * x1 - y1 (mod p)
        qc <<= qu.neg_ctrl(f1, f2, f3) @ minus_x_to_y_mod_n(y1, y3, self.p)
        qc <<= qu.neg_ctrl(f1, f2, f3) @ x_mul_y_add_zc_mod_n_to_zc(
            lam, x1, y3 + carry, self.p
        )
        return qc

    def _xx_add_yy_to_zz0_normal(
        self, xx: QubitTarget, yy: QubitTarget, zz: QubitTarget, ancilla: QubitTarget
    ) -> IGate:
        """
        doubling, xx+yy=O 無視バージョン
        """
        qc = CompositeGate()
        x1, y1 = xx[: self.bit_num], xx[self.bit_num : self.bit_num * 2]
        x2, y2 = yy[: self.bit_num], yy[self.bit_num : self.bit_num * 2]
        x3, y3 = zz[: self.bit_num], zz[self.bit_num : self.bit_num * 2]
        f1 = ancilla[0:1]
        f2 = ancilla[1:2]
        carry = ancilla[2:3]
        lam = ancilla[3 : self.bit_num + 3]
        aa = ancilla[self.bit_num + 3 :]

        # f1=(x1==0 && y1==0) (x1,y1)=O
        qc <<= qu.neg_ctrl(*x1, *y1) @ g.x(f1)
        # f2=(x2==0 && y2==0) (x2,y2)=O
        qc <<= qu.neg_ctrl(*x2, *y2) @ g.x(f2)
        # ctrl(f1) x3=x2
        qc <<= qu.ctrl(f1, x2) @ g.x(x3)
        # ctrl(f1) y3=y2
        qc <<= qu.ctrl(f1, y2) @ g.x(y3)
        # ctrl(f2) x3=x1
        qc <<= qu.ctrl(f2, x1) @ g.x(x3)
        # ctrl(f2) y3=y1
        qc <<= qu.ctrl(f2, y1) @ g.x(y3)
        # negctrl(f1) x3=-x2
        qc <<= qu.neg_ctrl(f1) @ minus_x_to_y_mod_n(x2, x3, self.p)
        # x2-=x1 (mod p)
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(x1, x2 + carry, self.p)
        # y2-=y1 (mod p)
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(
            y1,
            y2 + carry,
            self.p,
        )
        # aa=x2^-1 (mod p)
        qc <<= x_pow_a_mod_n_to_z_with_divs(
            x2,
            self.p - 2,
            self.p,
            aa[: self.bit_num],
            aa[self.bit_num :] + carry,
            self.inv_list,
        )
        # lambda=y2*aa (mod p)
        qc <<= x_mul_y_mod_n_to_zc(y2, aa, lam + carry, self.p)
        # negctrl(f1, f2, f3) x3 += lambda^2 - x1 (mod p)
        qc <<= qu.neg_ctrl(f1, f2) @ x_pow_2_add_zc_to_zc(lam, x3 + carry, self.p)
        qc <<= (
            qu.neg_ctrl(f1, f2) @ qu.inv @ x_add_zc_mod_n_to_zc(x1, x3 + carry, self.p)
        )
        # x1 -= x3 (mod p)
        qc <<= qu.inv @ x_add_zc_mod_n_to_zc(x3, x1 + carry, self.p)
        # negctrl(f1, f2) y3 = lambda * x1 - y1 (mod p)
        qc <<= qu.neg_ctrl(f1, f2) @ minus_x_to_y_mod_n(y1, y3, self.p)
        qc <<= qu.neg_ctrl(f1, f2) @ x_mul_y_add_zc_mod_n_to_zc(
            lam, x1, y3 + carry, self.p
        )
        return qc

    def _xx_add_Q_to_zz0_no_dbl(
        self, xx: QubitTarget, Q: Tuple[int, int], zz: QubitTarget, ancilla: QubitTarget
    ) -> IGate:
        """
        doubling無視バージョン
        """
        qc = CompositeGate()
        x1, y1 = Q[0], Q[1]
        x2, y2 = xx[: self.bit_num], xx[self.bit_num : self.bit_num * 2]
        x3, y3 = zz[: self.bit_num], zz[self.bit_num : self.bit_num * 2]
        f2 = ancilla[0:1]
        f3 = ancilla[1:2]
        carry = ancilla[2:3]
        lam = ancilla[3 : self.bit_num + 3]
        aa = ancilla[self.bit_num + 3 :]

        # f2=(x2==0 && y2==0) (x2,y2)=O
        qc <<= qu.neg_ctrl(*x2, *y2) @ g.x(f2)
        # ctrl(f2) x3=x1
        qc <<= qu.ctrl(f2) @ a_to_x(x1, x3)
        # ctrl(f2) y3=y1
        qc <<= qu.ctrl(f2) @ a_to_x(y1, y3)
        # x3=-x2
        qc <<= minus_x_to_y_mod_n(x2, x3, self.p)
        # x2-=x1 (mod p)
        qc <<= qu.inv @ xc_add_a_mod_n_to_xc(x2 + carry, x1, self.p)
        # y2-=y1 (mod p)
        qc <<= qu.inv @ xc_add_a_mod_n_to_xc(
            y2 + carry,
            y1,
            self.p,
        )
        # f3=(x2==0) 結果 O doublingなし
        qc <<= qu.neg_ctrl(*x2) @ g.x(f3)
        # aa=x2^-1 (mod p)
        qc <<= x_pow_a_mod_n_to_z_with_divs(
            x2,
            self.p - 2,
            self.p,
            aa[: self.bit_num],
            aa[self.bit_num :] + carry,
            self.inv_list,
        )
        # lambda=y2*aa (mod p)
        qc <<= x_mul_y_mod_n_to_zc(y2, aa, lam + carry, self.p)

        # negctrl(f2, f3) x3 += lambda^2 - x1 (mod p)
        qc <<= qu.neg_ctrl(f2, f3) @ x_pow_2_add_zc_to_zc(lam, x3 + carry, self.p)
        qc <<= qu.neg_ctrl(f2, f3) @ xc_add_a_mod_n_to_xc(
            x3 + carry, self.p - x1, self.p
        )
        # negctrl(f2, f3) y3 = lambda * x1 - lambda * x3 - y1 (mod p)
        qc <<= qu.neg_ctrl(f2, f3) @ x_mul_a_mod_n_to_zc(
            lam,
            x1,
            y3 + carry,
            self.p,
        )
        qc <<= (
            qu.neg_ctrl(f2, f3)
            @ qu.inv
            @ x_mul_y_add_zc_mod_n_to_zc(lam, x3, y3 + carry, self.p)
        )
        qc <<= qu.neg_ctrl(f2, f3) @ xc_add_a_mod_n_to_xc(
            y3 + carry, self.p - y1, self.p
        )
        return qc

    def xx_add_yy_to_zz(
        self, xx: QubitTarget, yy: QubitTarget, zz: QubitTarget, ancilla: QubitTarget
    ) -> IGate:
        if len(xx) < self.bit_num * 2:
            raise ValueError("x is too small")
        qc = CompositeGate()
        # 最後の計算用
        last_ancilla = ancilla[: self.bit_num + 1]
        # フラグ 3bit
        flag = ancilla[self.bit_num + 1 : self.bit_num + 4]
        # lambda
        lam = ancilla[self.bit_num + 4 : self.bit_num * 2 + 4]
        # ext
        ext = ancilla[self.bit_num * 2 + 4 :]

        qc <<= self.calc_xx_add_yy_to_flag_lambda(xx, yy, flag, lam, ext)

        # 最終座標
        qc <<= qu.neg_ctrl(*flag) @ self.xx_add_yy_with_lambda_to_zz(
            xx, yy, lam, qu.get_qubits(zz), last_ancilla
        )

        # xx=Oの処理
        qc <<= (
            qu.ctrl(flag[0])
            @ qu.ctrl(yy[: self.bit_num * 2])
            @ g.x(zz[: self.bit_num * 2])
        )
        # yy=Oの処理
        qc <<= (
            qu.ctrl(flag[1])
            @ qu.ctrl(xx[: self.bit_num * 2])
            @ g.x(zz[: self.bit_num * 2])
        )

        # 戻し
        qc <<= qu.inv @ self.calc_xx_add_yy_to_flag_lambda(xx, yy, flag, lam, ext)

        return qc

    def xx_add_Q_to_xx(
        self, xx: QubitTarget, Q: Tuple[int, int], ancilla: QubitTarget
    ) -> IGate:
        if len(xx) < self.bit_num * 2:
            raise ValueError("x is too small")
        qc = CompositeGate()
        qc <<= self.xx_add_Q_to_zz(
            xx, Q, ancilla[: self.bit_num * 2], ancilla[self.bit_num * 2 :]
        )
        qc <<= g.swap(xx, ancilla[: self.bit_num * 2])
        qc <<= qu.inv @ self.xx_add_Q_to_zz(
            xx,
            (Q[0], self.p - Q[1]),
            ancilla[: self.bit_num * 2],
            ancilla[self.bit_num * 2 :],
        )
        return qc

    def xx_add_yy_to_yy(
        self, xx: QubitTarget, yy: QubitTarget, ancilla: QubitTarget
    ) -> IGate:
        if len(xx) < self.bit_num * 2:
            raise ValueError("x is too small")
        qc = CompositeGate()
        qc <<= self.xx_add_yy_to_zz(
            xx, yy, ancilla[: self.bit_num * 2], ancilla[self.bit_num * 2 :]
        )
        qc <<= g.swap(yy, ancilla[: self.bit_num * 2])
        # xx.yをマイナスにする
        qc <<= minus_xc_mod_n_to_xc(
            xx[self.bit_num : self.bit_num * 2]
            + ancilla[self.bit_num * 2 : self.bit_num * 2 + 1],
            self.p,
        )
        # xx.y -=n, ~xx.y + 1
        qc <<= qu.inv @ self.xx_add_yy_to_zz(
            xx, yy, ancilla[: self.bit_num * 2], ancilla[self.bit_num * 2 :]
        )

        # -xx.y の戻し
        qc <<= qu.inv @ minus_xc_mod_n_to_xc(
            xx[self.bit_num : self.bit_num * 2] + ancilla[:1], self.p
        )

        return qc

    def x_mul_P_add_yy_to_yy(
        self, x: QubitTarget, P: Tuple[int, int], yy: QubitTarget, ancilla: QubitTarget
    ) -> IGate:
        qc = CompositeGate()
        pp = P
        for i in range(len(x)):
            print(f"2^{i} * {P} = {pp}")
            qc <<= qu.ctrl(x[i]) @ self.xx_add_Q_to_xx(yy, pp, ancilla)
            pp = self.add(pp, pp)
        return qc

    def _set_xy_PQ_to_zz(
        self, xy: QubitTarget, P: Tuple[int, int], Q: Tuple[int, int], zz: QubitTarget
    ) -> IGate:
        """
        xy=2bit
        """
        qc = CompositeGate()
        R = self.add(P, Q)
        # 0000
        # 0010 |01>
        # 0100 |10>
        # 0110 |01> |10>
        # 1000 |11>
        # 1010 |11> |01>
        # 1100 |11> |10>
        # 1110 |11> |01> |10>
        for i in range(self.bit_num * 2):
            iy, ix = divmod(i, self.bit_num)
            bit = 1 << ix
            qb = zz[i]
            flag = 0
            if P[iy] & bit:  # |01>
                flag |= 1
            if Q[iy] & bit:  # |10>
                flag |= 2
            if R[iy] & bit:  # |11>
                flag |= 4
            if flag == 3:
                qc <<= qu.ctrl(xy[0]) @ g.x(qb)
                qc <<= qu.ctrl(xy[1]) @ g.x(qb)
            elif flag == 5:
                qc <<= qu.ctrl(xy[0]) @ g.x(qb)
            elif flag == 6:
                qc <<= qu.ctrl(xy[1]) @ g.x(qb)
            elif flag == 7:
                qc <<= g.x(qb)
                qc <<= qu.neg_ctrl(*xy) @ g.x(qb)
            else:
                for j in range(3):
                    if flag & (1 << j):
                        if j == 0:
                            qc <<= qu.neg_ctrl(xy[1]) @ qu.ctrl(xy[0]) @ g.x(qb)
                        elif j == 1:
                            qc <<= qu.neg_ctrl(xy[0]) @ qu.ctrl(xy[1]) @ g.x(qb)
                        elif j == 2:
                            qc <<= qu.ctrl(*xy) @ g.x(qb)
        return qc

    def _set_xyz_PQR_to_zz(
        self,
        xyz: QubitTarget,
        P: Tuple[int, int],
        Q: Tuple[int, int],
        R: Tuple[int, int],
        zz: QubitTarget,
    ) -> IGate:
        """
        xyz=3bit
        """
        qc = CompositeGate()
        S = []
        for i in range(8):
            A = (0, 0)
            if i & 1:
                A = self.add(A, P)
            if i & 2:
                A = self.add(A, Q)
            if i & 4:
                A = self.add(A, R)
            S.append(A)
        # 0000
        # 0010 |01>
        # 0100 |10>
        # 0110 |01> |10>
        # 1000 |11>
        # 1010 |11> |01>
        # 1100 |11> |10>
        # 1110 |11> |01> |10>
        for i in range(self.bit_num * 2):
            iy, ix = divmod(i, self.bit_num)
            bit = 1 << ix
            qb = zz[i]
            one_flags = []
            for j in range(8):
                if S[j][iy] & bit:
                    one_flags.append(j)
            if len(one_flags) >= 6:
                # 反転した上で逆を使う
                one_flags = [j for j in range(8) if j not in one_flags]
                qc <<= g.x(qb)
            elif len(one_flags) > 2:
                # 各ビットに分解してチェック
                # 制御ビットが1つをコスト5として、掛け算していく
                # neg_ctrlは１つにつき 2 を足す
                # 単なる x は 1 を足す
                target = [0, 0, 0, [], []]
                for j in range(3):
                    zero_cnt = []
                    one_cnt = []
                    jbit = 1 << j
                    for flag in one_flags:
                        if flag & jbit:
                            one_cnt.append(flag)
                        else:
                            zero_cnt.append(flag)
                    cost = [0, 0]
                    for k, fl in enumerate([zero_cnt, one_cnt]):
                        if len(fl) == 1:
                            cost[k] = 3
                        elif len(fl) == 2:
                            if (fl[0] ^ fl[1]) | jbit == 7:
                                cost[k] = 6
                            else:
                                cost[k] = 2
                        elif len(fl) == 3:
                            cost[k] = 4
                        elif len(fl) == 4:
                            cost[k] = 1
                    if (
                        j == 0
                        or target[0] + target[1] > cost[0] + cost[1]
                        or (
                            target[0] + target[1] == cost[0] + cost[1]
                            and target[0] > cost[0]
                        )
                    ):
                        target = [cost[0], cost[1], j, zero_cnt, one_cnt]
                one_flags = []
                # ここで処理する
                t = target[2]
                tbit = 1 << t
                for k, flag in enumerate([target[3], target[4]]):
                    if len(flag) == 1:
                        one_flags.append(flag[0])
                    elif len(flag) == 2:
                        mask = flag[0] ^ flag[1]
                        if mask | tbit == 7:
                            # 全ビット制御
                            one_flags.append(flag[0])
                            one_flags.append(flag[1])
                        elif mask == 1:
                            qc <<= (flag[0] & 2 == 0 and qu.neg_ctrl(xyz[1]) or qu.ctrl(xyz[1])) @ (
                                flag[0] & 4 == 0 and qu.neg_ctrl(xyz[2]) or qu.ctrl(xyz[2])
                            ) @ g.x(qb)
                        elif mask == 2:
                            qc <<= (flag[0] & 1 == 0 and qu.neg_ctrl(xyz[0]) or qu.ctrl(xyz[0])) @ (
                                flag[0] & 4 == 0 and qu.neg_ctrl(xyz[2]) or qu.ctrl(xyz[2])
                            ) @ g.x(qb)
                        else:
                            qc <<= (flag[0] & 1 == 0 and qu.neg_ctrl(xyz[0]) or qu.ctrl(xyz[0])) @ (
                                flag[0] & 2 == 0 and qu.neg_ctrl(xyz[1]) or qu.ctrl(xyz[1])
                            ) @ g.x(qb)
                    elif len(flag) == 3:
                        yflg = tbit
                        if k == 0:
                            qc <<= qu.neg_ctrl(xyz[t]) @ g.x(qb)
                            yflg = 0
                        else:
                            qc <<= qu.ctrl(xyz[t]) @ g.x(qb)
                        for j in range(8):
                            if (j & tbit) != yflg:
                                continue
                            if j not in flag:
                                one_flags.append(j)
                                break
                    elif len(flag) == 4:
                        if k == 0:
                            qc <<= qu.neg_ctrl(xyz[t]) @ g.x(qb)
                        else:
                            qc <<= qu.ctrl(xyz[t]) @ g.x(qb)
            for flag in one_flags:
                qc <<= ctrl_with_state(xyz, flag) @ g.x(qb)
        return qc

    def x_mul_P_add_y_mul_Q_to_zz(
        self,
        x: QubitTarget,
        P: Tuple[int, int],
        y: QubitTarget,
        Q: Tuple[int, int],
        zz: QubitTarget,
        ancilla: QubitTarget,
    ) -> IGate:
        qc = CompositeGate()
        qc <<= self._set_xy_PQ_to_zz(x[0:1] + y[0:1], P, Q, zz)
        tmpqc = CompositeGate()
        for i in range(1, len(x)):
            P = self.add(P, P)
            Q = self.add(Q, Q)
            # print(i, P, Q)
            tmpqc <<= self._set_xy_PQ_to_zz(
                x[i : i + 1] + y[i : i + 1],
                P,
                Q,
                ancilla[(i - 1) * self.bit_num * 2 : i * self.bit_num * 2],
            )
        qc <<= tmpqc
        ext = ancilla[(len(x) - 1) * self.bit_num * 2 :]
        for i in range(1, len(x)):
            print(f"({i}/{len(x)})")
            qc <<= self.xx_add_yy_to_yy(
                ancilla[(i - 1) * self.bit_num * 2 : i * self.bit_num * 2], zz, ext
            )
        qc <<= qu.inv @ tmpqc
        print(f"({len(x)}/{len(x)})")
        return qc

    def x_mul_P_add_y_mul_Q_to_zz_v2(
        self,
        x: QubitTarget,
        P: Tuple[int, int],
        y: QubitTarget,
        Q: Tuple[int, int],
        zz: QubitTarget,
        ancilla: QubitTarget,
    ) -> IGate:
        qc = CompositeGate()
        tmpqc = CompositeGate()
        regs = []
        ix = 0
        for i in range(0, len(x)):
            # print(i, P, Q)
            reg = ancilla[ix : ix + self.bit_num * 2]
            regs.append(reg)
            ix += self.bit_num * 2
            tmpqc <<= self._set_xy_PQ_to_zz(
                x[i : i + 1] + y[i : i + 1],
                P,
                Q,
                reg,
            )
            P = self.add(P, P)
            Q = self.add(Q, Q)
        # 掛け算を続ける
        while len(regs) > 1:
            xx = regs.pop(0)
            print(f"Rest: {len(regs)}")
            yy = regs.pop(0)
            z = ancilla[ix : ix + self.bit_num * 2]
            ix += self.bit_num * 2
            anc = ancilla[ix : ix + self.get_add_ancilla_size()]
            ix += self.get_add_ancilla_size()
            regs.append(z)
            tmpqc <<= self._xx_add_yy_to_zz0(xx, yy, z, anc)
        qc <<= tmpqc
        qc <<= qu.ctrl(regs[0]) @ g.x(zz)
        print("Make Uncompute")
        qc <<= qu.inv @ tmpqc
        # qc <<= qu.ResetGate(ancilla)
        print("Complete aG+bQ")
        return qc

    def x_mul_P_add_y_mul_Q_to_zz_v4_wide(
        self,
        x: QubitTarget,
        P: Tuple[int, int],
        y: QubitTarget,
        Q: Tuple[int, int],
        zz: QubitTarget,
        ancilla: QubitTarget,
    ) -> IGate:
        qc = CompositeGate()
        tmpqc = CompositeGate()
        regs_P = []
        regs_Q = []
        ix = 0
        # ビット数によって分ける
        # 3の倍数: P,Qそれぞれ３組ペアずつ
        # 3の倍数+1: P,Qそれぞれ3組ペア、最上位ビットは P+Q の2つ
        # 3の倍数+2: P,Qそれぞれ3組ペアで、上位ビット2つでそれぞれペア
        nP, nQ = P, Q
        st = 0
        if len(x) % 3 == 2:
            # あまりは最上位と最下位で足す
            st = 1
            nP = self.add(nP, nP)
            nQ = self.add(nQ, nQ)
        for a in range(len(x) // 3):
            i = st + a * 3
            # print(i, P, Q)
            plist = []
            qlist = []
            for _ in range(3):
                plist.append(nP)
                qlist.append(nQ)
                nP = self.add(nP, nP)
                nQ = self.add(nQ, nQ)
            reg = ancilla[ix : ix + self.bit_num * 2]
            regs_P.append(reg)
            ix += self.bit_num * 2
            print("set P", plist)
            tmpqc <<= self._set_xyz_PQR_to_zz(
                x[i : i + 3], plist[0], plist[1], plist[2], reg
            )
            reg = ancilla[ix : ix + self.bit_num * 2]
            regs_Q.append(reg)
            ix += self.bit_num * 2
            print("set Q", qlist)
            tmpqc <<= self._set_xyz_PQR_to_zz(
                y[i : i + 3],
                qlist[0],
                qlist[1],
                qlist[2],
                reg,
            )
        ext_P = None
        ext_Q = None
        ext_PQ = None
        # あまりで処理
        # 0: regs_Pの最後をext_Pに移動
        if len(x) % 3 == 0:
            ext_P = regs_P.pop(-1)
            ext_Q = regs_Q.pop(-1)
        elif len(x) % 3 == 1:
            ext_PQ = ancilla[ix : ix + self.bit_num * 2]
            ix += self.bit_num * 2
            tmpqc <<= self._set_xy_PQ_to_zz(
                x[-1:] + y[-1:],
                nP,
                nQ,
                ext_PQ,
            )
        elif len(x) % 3 == 2:
            ext_P = ancilla[ix : ix + self.bit_num * 2]
            ix += self.bit_num * 2
            tmpqc <<= self._set_xy_PQ_to_zz(
                x[:1] + x[-1:],
                P,
                nP,
                ext_P,
            )
            ext_Q = ancilla[ix : ix + self.bit_num * 2]
            ix += self.bit_num * 2
            tmpqc <<= self._set_xy_PQ_to_zz(
                y[:1] + y[-1:],
                Q,
                nQ,
                ext_Q,
            )
        last_reg = ancilla[ix : ix + self.bit_num * 2]
        ix += self.bit_num * 2

        # 掛け算を続ける
        while len(regs_P) > 1:
            xx = regs_P.pop(0)
            print(f"add P(normal): {len(regs_P)}")
            yy = regs_P.pop(0)
            z = ancilla[ix : ix + self.bit_num * 2]
            ix += self.bit_num * 2
            anc = ancilla[ix : ix + self.get_add_ancilla_size()]
            ix += self.get_add_ancilla_size()
            regs_P.append(z)
            tmpqc <<= self._xx_add_yy_to_zz0_normal(xx, yy, z, anc)
        while len(regs_Q) > 1:
            xx = regs_Q.pop(0)
            print(f"add Q(normal): {len(regs_Q)}")
            yy = regs_Q.pop(0)
            z = ancilla[ix : ix + self.bit_num * 2]
            ix += self.bit_num * 2
            anc = ancilla[ix : ix + self.get_add_ancilla_size()]
            ix += self.get_add_ancilla_size()
            regs_Q.append(z)
            tmpqc <<= self._xx_add_yy_to_zz0_normal(xx, yy, z, anc)

        if ext_P is not None and ext_Q is not None:
            if len(regs_P) > 0:
                print("add P(last no double),add Q(last no double)")
                zP = ancilla[ix : ix + self.bit_num * 2]
                ix += self.bit_num * 2
                anc = ancilla[ix : ix + self.get_add_ancilla_size()]
                ix += self.get_add_ancilla_size()
                tmpqc <<= self._xx_add_yy_to_zz0_no_dbl(regs_P[0], ext_P, zP, anc)
                zQ = ancilla[ix : ix + self.bit_num * 2]
                ix += self.bit_num * 2
                anc = ancilla[ix : ix + self.get_add_ancilla_size()]
                ix += self.get_add_ancilla_size()
                tmpqc <<= self._xx_add_yy_to_zz0_no_dbl(regs_Q[0], ext_Q, zQ, anc)
                # 最後に足し算する
                tmpqc <<= self._xx_add_yy_to_zz0_no_dbl(
                    zP, zQ, last_reg, ancilla[ix : ix + self.get_add_ancilla_size()]
                )
            else:
                print("add P(all)+Q(all)")
                tmpqc <<= self._xx_add_yy_to_zz0(
                    ext_P,
                    ext_Q,
                    last_reg,
                    ancilla[ix : ix + self.get_add_ancilla_size()],
                )
        elif ext_PQ is not None:
            print("add P(last)+Q'last)")
            zPQ = ancilla[ix : ix + self.bit_num * 2]
            ix += self.bit_num * 2
            anc = ancilla[ix : ix + self.get_add_ancilla_size()]
            ix += self.get_add_ancilla_size()
            tmpqc <<= self._xx_add_yy_to_zz0(regs_P[0], regs_Q[0], zPQ, anc)
            # 最後に足し算する
            tmpqc <<= self._xx_add_yy_to_zz0(
                zPQ, ext_PQ, last_reg, ancilla[ix : ix + self.get_add_ancilla_size()]
            )

        qc <<= tmpqc
        qc <<= qu.ctrl(last_reg) @ g.x(zz)
        print("Make Uncompute")
        qc <<= qu.inv @ tmpqc
        # qc <<= qu.ResetGate(ancilla)
        print("Complete aG+bQ")
        return qc

    def x_mul_P_add_y_mul_Q_to_zz_v3(
        self,
        x: QubitTarget,
        P: Tuple[int, int],
        y: QubitTarget,
        Q: Tuple[int, int],
        zz: QubitTarget,
        ancilla: QubitTarget,
    ) -> IGate:
        qc = CompositeGate()
        org_P, org_Q = P, Q
        qc <<= self._set_xy_PQ_to_zz(x[0:1] + y[0:1], P, Q, zz)
        tmp_qb = ancilla[: self.bit_num * 2]
        nxt_qb = ancilla[self.bit_num * 2 : self.bit_num * 4]
        tmp_zz = ancilla[self.bit_num * 4 : self.bit_num * 6]
        anc = ancilla[self.bit_num * 6 :]
        for i in range(1, len(x)):
            print(f"({i}/{len(x)})")
            P = self.add(P, P)
            Q = self.add(Q, Q)
            # print(i, P, Q)
            qc <<= self._set_xy_PQ_to_zz(
                x[i : i + 1] + y[i : i + 1],
                P,
                Q,
                tmp_qb,
            )
            qc <<= self._xx_add_yy_to_zz0(zz, tmp_qb, tmp_zz, anc)
            qc <<= qu.ctrl(tmp_zz) @ g.x(nxt_qb)
            qc <<= qu.inv @ self._xx_add_yy_to_zz0(zz, tmp_qb, tmp_zz, anc)
            qc <<= g.swap(nxt_qb, zz)
            if i == 1:
                qc <<= qu.inv @ self._set_xy_PQ_to_zz(
                    x[0:1] + y[0:1], org_P, org_Q, nxt_qb
                )
            else:
                # 引き算
                qc <<= minus_xc_mod_n_to_xc(
                    tmp_qb[self.bit_num : self.bit_num * 2] + anc,
                    self.p,
                )
                qc <<= self._xx_add_yy_to_zz0(zz, tmp_qb, tmp_zz, anc)
                qc <<= qu.ctrl(tmp_zz) @ g.x(nxt_qb)
                qc <<= qu.inv @ self._xx_add_yy_to_zz0(zz, tmp_qb, tmp_zz, anc)
                qc <<= qu.inv @ minus_xc_mod_n_to_xc(
                    tmp_qb[self.bit_num : self.bit_num * 2] + anc,
                    self.p,
                )
            qc <<= qu.inv @ self._set_xy_PQ_to_zz(
                x[i : i + 1] + y[i : i + 1],
                P,
                Q,
                tmp_qb,
            )

        # qc <<= qu.ResetGate(ancilla)
        print("Complete aG+bQ")
        return qc

    def x_mul_P_add_y_mul_Q_to_zz_v4_compact(
        self,
        x: QubitTarget,
        P: Tuple[int, int],
        y: QubitTarget,
        Q: Tuple[int, int],
        zz: QubitTarget,
        ancilla: QubitTarget,
    ) -> IGate:
        qc = CompositeGate()
        pos_list = []
        nP = P
        xy = x[:] + y[:]
        for i in range(len(x) * 2):
            if i == len(x):
                nP = Q
            pos_list.append(nP)
            nP = self.add(nP, nP)
        tmp_qb = ancilla[: self.bit_num * 2]
        nxt_qb = ancilla[self.bit_num * 2 : self.bit_num * 4]
        tmp_zz = ancilla[self.bit_num * 4 : self.bit_num * 6]
        anc = ancilla[self.bit_num * 6 :]
        # pos_listを3つのグループに分ける
        group_list = []
        for i in range(math.ceil(len(pos_list) / 3)):
            group_list.append(pos_list[i * 3 : i * 3 + 3])
        if len(pos_list) % 6 == 0 and len(group_list) > 2:
            # Pの合成について先頭グループと最後のグループを調整する
            last_ix = len(group_list) // 2 - 1
            group_list[0][2], group_list[last_ix][2] = (
                group_list[last_ix][2],
                group_list[0][2],
            )
            xy[2], xy[last_ix * 3 + 2] = xy[last_ix * 3 + 2], xy[2]
        elif len(x) % 3 == 2:
            # xをnodblで実行できる
            last_ix = len(group_list) // 2 - 1
            qpos = group_list[last_ix].pop(-1)
            group_list[-1].append(qpos)
            xy = (
                xy[: last_ix * 3 + 2]
                + xy[last_ix * 3 + 3 :]
                + xy[last_ix * 3 + 2 : last_ix * 3 + 3]
            )
            group_list[0][0], group_list[last_ix][0] = (
                group_list[last_ix][0],
                group_list[0][0],
            )
            xy[0], xy[last_ix * 3] = xy[last_ix * 3], xy[0]
        ix = 0
        for group in group_list:
            if ix == 0:
                print("set", group)
                qc <<= self._set_xyz_PQR_to_zz(
                    xy[0:3], group[0], group[1], group[2], zz
                )
                ix = 3
                continue
            else:
                sz = len(group)
                label = ""
                add_proc = self._xx_add_yy_to_zz0
                sub_proc = self._xx_add_yy_to_zz0
                if sz == 1:
                    label = "const"
                elif ix + len(group) < len(x):
                    # 最後でなければ通常足し算のみ
                    add_proc = self._xx_add_yy_to_zz0_normal
                    sub_proc = self._xx_add_yy_to_zz0_no_dbl
                    label = "normal"
                elif ix + len(group) == len(x):
                    # doublingなし
                    add_proc = self._xx_add_yy_to_zz0_no_dbl
                    label = "no double"
                print("add", "%d/%d" % (ix, len(pos_list)), group, label)
                # 座標の数に応じて処理を変更
                if sz == 1:
                    # 定数座標の足し算
                    qc <<= qu.ctrl(xy[ix]) @ self.xx_add_Q_to_zz0(
                        zz, group[0], tmp_zz, anc
                    )
                    qc <<= qu.ctrl(xy[ix]) @ qu.ctrl(tmp_zz) @ g.x(nxt_qb)
                    qc <<= (
                        qu.ctrl(xy[ix])
                        @ qu.inv
                        @ self.xx_add_Q_to_zz0(zz, group[0], tmp_zz, anc)
                    )
                    qc <<= qu.ctrl(xy[ix]) @ g.swap(nxt_qb, zz)
                    qc <<= qu.ctrl(xy[ix]) @ self.xx_add_Q_to_zz0(
                        zz, (group[0][0], self.p - group[0][1]), tmp_zz, anc
                    )
                    qc <<= qu.ctrl(xy[ix]) @ qu.ctrl(tmp_zz) @ g.x(nxt_qb)
                    qc <<= (
                        qu.ctrl(xy[ix])
                        @ qu.inv
                        @ self.xx_add_Q_to_zz0(
                            zz, (group[0][0], self.p - group[0][1]), tmp_zz, anc
                        )
                    )
                    ix += 1
                    continue
                elif sz == 2:
                    # 2つの座標の合成
                    qc <<= self._set_xy_PQ_to_zz(
                        xy[ix : ix + 2],
                        group[0],
                        group[1],
                        tmp_qb,
                    )
                elif sz == 3:
                    # 3つの座標の合成
                    qc <<= self._set_xyz_PQR_to_zz(
                        xy[ix : ix + 3],
                        group[0],
                        group[1],
                        group[2],
                        tmp_qb,
                    )

                # 足し算する
                qc <<= add_proc(zz, tmp_qb, tmp_zz, anc)
                qc <<= qu.ctrl(tmp_zz) @ g.x(nxt_qb)
                qc <<= qu.inv @ add_proc(zz, tmp_qb, tmp_zz, anc)
                qc <<= g.swap(nxt_qb, zz)
                qc <<= minus_xc_mod_n_to_xc(
                    tmp_qb[self.bit_num : self.bit_num * 2] + anc,
                    self.p,
                )
                qc <<= sub_proc(zz, tmp_qb, tmp_zz, anc)
                qc <<= qu.ctrl(tmp_zz) @ g.x(nxt_qb)
                qc <<= qu.inv @ sub_proc(zz, tmp_qb, tmp_zz, anc)
                qc <<= qu.inv @ minus_xc_mod_n_to_xc(
                    tmp_qb[self.bit_num : self.bit_num * 2] + anc,
                    self.p,
                )
                # 戻し
                if sz == 2:
                    qc <<= qu.inv @ self._set_xy_PQ_to_zz(
                        xy[ix : ix + 2],
                        group[0],
                        group[1],
                        tmp_qb,
                    )
                    ix += 2
                elif sz == 3:
                    qc <<= qu.inv @ self._set_xyz_PQR_to_zz(
                        xy[ix : ix + 3],
                        group[0],
                        group[1],
                        group[2],
                        tmp_qb,
                    )
                    ix += 3
        # qc <<= qu.ResetGate(ancilla)
        print("Complete aG+bQ")
        return qc


if __name__ == "__main__":
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    Gx = int(sys.argv[3])
    Gy = int(sys.argv[4])
    p = int(sys.argv[5])
    ecc = EccQuantum(a, b, (Gx, Gy), p)
    print(f"a={a}, b={b}, G=({Gx},{Gy}), p={p}, order={ecc.order}")
    n = ecc.order
    for i in range(n + 1):
        P = ecc.get_nG(i)
        print(f"{i}G = {P}, valid={ecc.is_valid(P)}")
