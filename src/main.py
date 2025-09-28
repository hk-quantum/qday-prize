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

import numpy as np
from qiskit import QuantumCircuit, transpile, QuantumRegister, ClassicalRegister
from typing import Tuple, Dict
import hkutils.quantum as qu
from hkutils.ecc import EccQuantum
import hkutils.stdgates as g
import hkutils.backend as sim_mod
import hkutils.ecc as ecc_mod
import math
import time
import ibm_backend
import json
import sys


def get_ecc(bit_num: int) -> Tuple[EccQuantum, Tuple[int, int]]:
    with open("data/curves.json", "r", encoding="utf-8") as f:
        data = json.load(f)  # JSONをPythonオブジェクトに変換
    for dt in data:
        if dt["bit_length"] == bit_num:
            ecc = EccQuantum(0, 7, tuple(dt["generator_point"]), dt["prime"])
            return ecc, tuple(dt["public_key"])
    raise ValueError(f"Invalid Bit Size {bit_num}")


def check_counts(counts: Dict[str, int], ecc: EccQuantum, Q: Tuple[int, int]):
    def is_ok(x, y) -> int:
        try:
            d = (pow(x, -1, ecc.order)) * y % ecc.order
            if ecc.get_nG(d) == Q:
                return d
        except:
            pass
        return 0

    ok_count = 0
    result_d = 0
    sorted_items = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    for key, cnt in sorted_items:
        vals = key.split()
        a = int(vals[0], 2)
        b = int(vals[1], 2)
        x = a * ecc.order / (1 << ecc.order.bit_length())
        y = b * ecc.order / (1 << ecc.order.bit_length())
        if (d := is_ok(math.floor(x), math.floor(y))) > 0:
            ok_count += cnt
        elif (d := is_ok(math.ceil(x), math.floor(y))) > 0:
            ok_count += cnt
        elif (d := is_ok(math.floor(x), math.ceil(y))) > 0:
            ok_count += cnt
        elif (d := is_ok(math.ceil(x), math.ceil(y))) > 0:
            ok_count += cnt
        else:
            d = 0
        if d > 0:
            print(f"*OK: <{key}>={cnt}, d={d}")
            result_d = d
        else:
            print(f"-NG: <{key}>={cnt}")
    print(f"Success: d={result_d} count={ok_count}")


if __name__ == "__main__":
    np.random.seed(None)
    bit_num = 4
    if len(sys.argv) > 1:
        bit_num = int(sys.argv[1])

    # https://www.qdayprize.org/curves.txt
    ecc, Q = get_ecc(bit_num)

    print(f"bit_length={bit_num} prime={ecc.p} G={ecc.G} n={ecc.order} Q={Q}")

    start_time = time.perf_counter()
    print("Creating Quantum Circuit")
    a_qb = QuantumRegister(ecc.order.bit_length(), name="qa")
    b_qb = QuantumRegister(ecc.order.bit_length(), name="qb")
    x_qb = QuantumRegister(ecc.bit_num, name="qx")
    y_qb = QuantumRegister(ecc.bit_num, name="qy")
    ancilla = QuantumRegister(
        ecc.bit_num * ((ecc.order - 2).bit_length() + 7) + 4, name="ancilla"
    )
    a_c = ClassicalRegister(len(a_qb), name="a")
    b_c = ClassicalRegister(len(b_qb), name="b")
    x_c = ClassicalRegister(len(x_qb), name="x")
    y_c = ClassicalRegister(len(y_qb), name="y")
    circuit = QuantumCircuit(a_qb, b_qb, x_qb, y_qb, ancilla, y_c, x_c, b_c, a_c)
    qc = qu.QiskitWrapperCircuit(circuit)
    qc <<= g.h(a_qb)
    qc <<= g.h(b_qb)

    qc <<= ecc.x_mul_P_add_yy_to_yy(
        a_qb, ecc.G, qu.get_qubits(x_qb) + qu.get_qubits(y_qb), ancilla
    )
    qc <<= ecc.x_mul_P_add_yy_to_yy(
        b_qb, Q, qu.get_qubits(x_qb) + qu.get_qubits(y_qb), ancilla
    )
    qc.barrier()
    qc.measure(x_qb, x_c)
    qc.measure(y_qb, y_c)
    qc.barrier()
    qc <<= ecc_mod.qft_dagger(a_qb)
    qc <<= ecc_mod.qft_dagger(b_qb)
    qc.measure(a_qb, a_c)
    qc.measure(b_qb, b_c)
    tm = time.perf_counter() - start_time
    print(f"Complete Quantum Circuit {tm:.3f}[s]")
    print(f"  qubit_size={qc.circuit.num_qubits} gate_count={len(qc.circuit.data)}")

    # 実行する
    # backend, t_qc = ibm_backend.get_backend(circuit)
    backend, t_qc = sim_mod.get_backend(circuit)
    jobs = backend.run([t_qc], shots=100)
    results = jobs.result()
    counts = results.get_counts()
    # キー順に表示
    check_counts(counts, ecc, Q)  # type:ignore
