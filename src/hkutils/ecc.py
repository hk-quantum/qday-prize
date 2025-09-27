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

from typing import Tuple
from qiskit.circuit import Qubit
from . import stdgates as g
from . import quantum as qu
from .quantum import QubitTarget, IGate, CompositeGate
import numpy as np
import time

"""
x3
y3
carry
flag(xx=O, xx+Q=O)
lambdaEccQuantum
dy
dx
dx^-1 (+dx^2,dx^4...tmp)
lambda * x3
flag(xx=O)
"""

def qft_dagger(reg: QubitTarget) -> qu.IGate:
    qc = CompositeGate()
    n = len(reg)
    for ix in range(n // 2):
        qc <<= g.swap(reg[ix], reg[n - ix - 1])
    for j in range(n):
        for m in range(j):
            qc <<= g.cp(-np.pi / float(2 ** (j - m)))(reg[m], reg[j])
            #qc <<= qu.ctrl(reg[m]) @ g.p(-np.pi / float(2 ** (j - m)))(reg[j])
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
    qc <<= x_add_a_to_x(xc[: bit_num + 1], a)
    qc <<= qu.inv @ x_add_a_to_x(xc[: bit_num + 1], n)
    qc <<= qu.ctrl(xc[bit_num]) @ x_add_a_to_x(xc[:bit_num], n)
    qc <<= qu.inv @ x_add_a_to_x(xc[: bit_num + 1], a)
    qc <<= x_add_a_to_x(xc[:bit_num], a)
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


def x_mul_a_mod_n_to_zc(x: QubitTarget, a: int, zc: QubitTarget, n: int) -> IGate:
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
    xをyで掛け算してnで割った余りをzに格納する
    zはbit_num + 1
    """
    bit_num = n.bit_length()
    if len(zc) < bit_num + 1:
        raise ValueError("zのビット数が足りません")
    qc = qu.CompositeGate()
    for i in reversed(range(bit_num)):
        for j in reversed(range(bit_num)):
            qc <<= qu.ctrl(x[i], y[j]) @ xc_add_a_mod_n_to_xc(
                zc[: bit_num + 1], ((1 << i) * (1 << j)) % n, n
            )
    return qc


def x_pow_2_mod_n_to_zc(x: QubitTarget, n: int, zc: QubitTarget) -> IGate:
    """
    xを2乗して n で割ってzに格納する
    zはbit_num + 1
    """
    bit_num = n.bit_length()
    qc = qu.CompositeGate()
    if len(zc) < bit_num + 1:
        raise ValueError("zのビット数が足りません")
    tz = zc[: bit_num * 2]
    for i in reversed(range(bit_num)):
        for j in reversed(range(bit_num)):
            if i == j:
                qc <<= qu.ctrl(x[i]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], ((1 << i) * (1 << j)) % n, n
                )
            else:
                qc <<= qu.ctrl(x[i], x[j]) @ xc_add_a_mod_n_to_xc(
                    zc[: bit_num + 1], ((1 << i) * (1 << j)) % n, n
                )
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
        lambdaの分母 y2-y1
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

    def debug_init(self, x: QubitTarget, y: QubitTarget, ancilla: QubitTarget):
        qc = CompositeGate()
        # Initialize the debug circuit
        qc <<= g.h(ancilla[: self.bit_num])
        for i in range(1, self.order):
            pos = self.get_nG(i)
            print("init", i, pos)
            mod = None
            for j in reversed(range(self.bit_num)):
                if i & (1 << j):
                    if mod is None:
                        mod = qu.ctrl(ancilla[j])
                    else:
                        mod = mod @ qu.ctrl(ancilla[j])
                else:
                    if mod is None:
                        mod = qu.neg_ctrl(ancilla[j])
                    else:
                        mod = mod @ qu.neg_ctrl(ancilla[j])
            qc <<= mod @ x_add_a_to_x(x, pos[0])  # type: ignore
            qc <<= mod @ x_add_a_to_x(y, pos[1])  # type: ignore
        return qc
