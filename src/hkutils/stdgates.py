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

from qiskit.circuit.library import (
    U1Gate,
    XGate,
    YGate,
    ZGate,
    HGate,
    SwapGate,
    SGate,
    SdgGate,
    TGate,
    TdgGate,
    SXGate,
    RXGate,
    RYGate,
    RZGate,
    CXGate,
    CYGate,
    CZGate,
    CPhaseGate,
    CRXGate,
    CRYGate,
    CRZGate,
    CHGate,
    CU3Gate,
    CCXGate,
    CSwapGate,
    U2Gate,
    U3Gate,
)
from .quantum import GateDef


def p(phase: float) -> GateDef:
    return GateDef(U1Gate(phase))


x = GateDef(XGate())
y = GateDef(YGate())
z = GateDef(ZGate())
h = GateDef(HGate())
s = GateDef(SGate())
sdg = GateDef(SdgGate())
t = GateDef(TGate())
tdg = GateDef(TdgGate())
sx = GateDef(SXGate())


def rx(theta: float) -> GateDef:
    return GateDef(RXGate(theta))


def ry(theta: float) -> GateDef:
    return GateDef(RYGate(theta))


def rz(theta: float) -> GateDef:
    return GateDef(RZGate(theta))


cx = GateDef(CXGate())
cy = GateDef(CYGate())
cz = GateDef(CZGate())


def cp(theta: float) -> GateDef:
    return GateDef(CPhaseGate(theta))


def crx(theta: float) -> GateDef:
    return GateDef(CRXGate(theta))


def cry(theta: float) -> GateDef:
    return GateDef(CRYGate(theta))


def crz(theta: float) -> GateDef:
    return GateDef(CRZGate(theta))


ch = GateDef(CHGate())


def cu(theta: float, phi: float, lam: float, gammba: float) -> GateDef:
    # TODO gamma
    return GateDef(CU3Gate(theta, phi, lam))


swap = GateDef(SwapGate())

ccx = GateDef(CCXGate())
cswap = GateDef(CSwapGate())
CX = cx


def phase(lam: float) -> GateDef:
    return p(lam)


def cphase(lam: float) -> GateDef:
    return cp(lam)


def u1(lam: float) -> GateDef:
    return GateDef(U1Gate(lam))


def u2(phi: float, lam: float) -> GateDef:
    return GateDef(U2Gate(phi, lam))


def u3(theta: float, phi: float, lam: float) -> GateDef:
    return GateDef(U3Gate(theta, phi, lam))
