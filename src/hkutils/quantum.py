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

from typing import List, Tuple, Union, overload
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit import Gate, Qubit
from abc import ABC, abstractmethod

QubitTarget = Union[QuantumRegister, Qubit, List[Qubit]]


def get_qubits(target: QubitTarget) -> List[Qubit]:
    if isinstance(target, QuantumRegister):
        return target[:]
    elif isinstance(target, Qubit):
        return [target]
    return target


class GateModifier:
    _ctrl: List[List[Qubit]]
    _neg_ctrl: List[List[Qubit]]
    _power: int

    def __init__(
        self, ctrl: List[QubitTarget], neg_ctrl: List[QubitTarget], power: int = 1
    ):
        self._ctrl = [get_qubits(c) for c in ctrl]
        self._neg_ctrl = [get_qubits(c) for c in neg_ctrl]
        self._power = power

    def merge(self, other: "GateModifier") -> "GateModifier":
        merged_ctrl = self._ctrl + other._ctrl
        merged_neg_ctrl = self._neg_ctrl + other._neg_ctrl
        merged_power = self._power * other._power
        return GateModifier(merged_ctrl, merged_neg_ctrl, merged_power)

    @overload
    def __matmul__(self, other: "GateModifier") -> "GateModifier": ...  # 実装は書かない
    @overload
    def __matmul__(self, other: "IGate") -> "IGate": ...  # 実装は書かない
    def __matmul__(
        self, other: Union["GateModifier", "IGate"]
    ) -> Union["GateModifier", "IGate"]:
        if isinstance(other, IGate):
            return other.apply(self)
        return self.merge(other)


class IGate(ABC):
    @abstractmethod
    def apply(self, modifier: GateModifier) -> "IGate":
        pass

    @abstractmethod
    def extract(self, qc: QuantumCircuit) -> None:
        pass


class QuGate(IGate):
    _gate: "GateDef"
    _qubits: List[List[Qubit]]
    _modifier: Union[GateModifier, None]

    def __init__(self, gate: "GateDef", *qubits: List[QubitTarget]):
        self._gate = gate
        self._qubits = [get_qubits(qb) for qb in qubits]
        self._modifier = None

    def apply(self, modifier: GateModifier) -> "QuGate":
        result = QuGate(self._gate, *self._qubits)
        if self._modifier is not None:
            result._modifier = self._modifier.merge(modifier)
        else:
            result._modifier = modifier
        return result

    def extract(self, qc: QuantumCircuit) -> None:
        bit_sz = 1

        def check_bit_size(qubits: List[List[Qubit]]):
            nonlocal bit_sz
            for qb in qubits:
                if len(qb) != 1:
                    if bit_sz == 1:
                        bit_sz = len(qb)
                    elif bit_sz != len(qb):
                        raise ValueError("Inconsistent qubit sizes")

        check_bit_size(self._qubits)
        if self._modifier is not None:
            check_bit_size(self._modifier._ctrl)
            check_bit_size(self._modifier._neg_ctrl)
        for bx in range(bit_sz):
            qubits = [qb[bx % len(qb)] for qb in self._qubits]
            ctrl = []
            neg_ctrl = []
            gate = self._gate
            if self._modifier is not None:
                ctrl = [c[bx % len(c)] for c in self._modifier._ctrl]
                neg_ctrl = [nc[bx % len(nc)] for nc in self._modifier._neg_ctrl]
                if self._modifier._power < 0:
                    gate = gate.inverse()
                if abs(self._modifier._power) > 1:
                    gate = gate.power(abs(self._modifier._power))
            self._gate.extract(qc, qubits, ctrl, neg_ctrl)


class CompositeGate(IGate):
    _gates: List[IGate]

    def __init__(self, *gates: IGate):
        self._gates = list(gates)

    def apply(self, modifier: GateModifier) -> "CompositeGate":
        if modifier._power == 0:
            return CompositeGate()
        if (
            len(modifier._ctrl) == 0
            and len(modifier._neg_ctrl) == 0
            and modifier._power == 1
        ):
            base = self._gates
        else:
            base = [
                gate.apply(GateModifier(modifier._ctrl, modifier._neg_ctrl, 1))
                for gate in self._gates
            ]
        if modifier._power < 0:
            base = list(reversed([inv @ gate for gate in self._gates]))
        return CompositeGate(*(base * abs(modifier._power)))

    def extract(self, qc: QuantumCircuit) -> None:
        for gate in self._gates:
            gate.extract(qc)

    def __lshift__(self, other: IGate):
        self.append(other)
        return self

    def __ilshift__(self, other: IGate):
        self.append(other)
        return self

    def append(self, gate: IGate) -> "CompositeGate":
        self._gates.append(gate)
        return self


class GateDef:
    _gate: Gate

    def __init__(self, gate: Union[Gate, QuantumCircuit]):
        if isinstance(gate, QuantumCircuit):
            self._gate = gate.to_gate()
        else:
            self._gate = gate

    def extract(
        self,
        qc: QuantumCircuit,
        qubits: List[Qubit],
        ctrl: List[Qubit],
        neg_ctrl: List[Qubit],
    ):
        qb_index = [qc.qubits.index(qb) for qb in qubits]
        ctrl_index = [qc.qubits.index(qb) for qb in ctrl]
        neg_ctrl_index = [qc.qubits.index(qb) for qb in neg_ctrl]
        if len(ctrl_index) > 0 or len(neg_ctrl_index) > 0:
            qc.append(
                self._gate.control(
                    len(ctrl_index) + len(neg_ctrl_index),
                    ctrl_state=(1 << len(ctrl_index)) - 1,
                ),
                ctrl_index + neg_ctrl_index + qb_index,
            )
        else:
            qc.append(self._gate, qb_index)

    def power(self, exponent: int) -> "GateDef":
        return GateDef(self._gate.power(exponent))  # type: ignore

    def inverse(self) -> "GateDef":
        return GateDef(self._gate.inverse())  # type: ignore

    def apply(self, *args: List[QubitTarget]) -> IGate:
        return QuGate(self, *args)

    def __call__(self, *args: List[QubitTarget]) -> IGate:
        return self.apply(*args)


class QiskitWrapperCircuit:
    _circuit: QuantumCircuit

    def __init__(self, circuit: QuantumCircuit):
        self._circuit = circuit

    @property
    def circuit(self) -> QuantumCircuit:
        return self._circuit

    def append(self, gate: IGate) -> "QiskitWrapperCircuit":
        # 適用
        gate.extract(self._circuit)
        return self

    def __ilshift__(self, other: IGate) -> "QiskitWrapperCircuit":
        return self.append(other)

    def barrier(self, *qargs: QubitTarget, label=None) -> "QiskitWrapperCircuit":
        self._circuit.barrier(qargs, label=label)
        return self

    def measure(
        self, qbits: QubitTarget, cregs: ClassicalRegister
    ) -> "QiskitWrapperCircuit":
        self._circuit.measure(qbits, cregs)
        return self


inv = GateModifier([], [], -1)


def ctrl(*targets: QubitTarget) -> GateModifier:
    return GateModifier(list(targets), [], 1)


def neg_ctrl(*targets: QubitTarget) -> GateModifier:
    return GateModifier([], list(targets), 1)
