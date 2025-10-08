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
import sys


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


def minus_xc_mod_n_to_xc(xc: QubitTarget, n: int) -> IGate:
    qc = CompositeGate()
    bit_num = n.bit_length()
    if len(xc) < bit_num + 1:
        raise ValueError("x is too small")
    qc <<= g.x(xc)
    qc <<= x_inc_to_x(xc)
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


def cx_mul_2_mod_n_to_xc(xc: QubitTarget, n: int) -> IGate:
    qc = CompositeGate()
    bit_num = n.bit_length()
    if len(xc) < bit_num + 1:
        raise ValueError("x is too small")

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

    def is_valid(self, P: Tuple[int, int]) -> bool:
        if P == (0, 0):
            return True
        x, y = P
        return (y * y - (x * x * x + self.a * x + self.b)) % self.p == 0

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
            if P[iy] & bit:  # |01>
                qc <<= qu.neg_ctrl(xy[1]) @ qu.ctrl(xy[0]) @ g.x(qb)
            if Q[iy] & bit:  # |10>
                qc <<= qu.ctrl(xy[1]) @ qu.neg_ctrl(xy[0]) @ g.x(qb)
            if R[iy] & bit:
                qc <<= qu.ctrl(*xy) @ g.x(qb)
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
