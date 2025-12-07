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

from qiskit.providers import BackendV2, Options, JobV1, JobStatus
from qiskit.result import Result
import numpy as np

from qiskit.circuit import ControlledGate
from qiskit.circuit.library import Measure, Reset, Barrier
from qiskit.result.models import ExperimentResultData, ExperimentResult

from qiskit import QuantumCircuit
import numpy as np
import uuid
from qiskit.circuit.library import Measure
from typing import List, Tuple, Any, Union, Dict, Generator
from numpy.typing import NDArray
import uuid
from .simulator import QuantumSimulator as CpuQuantumSimulator
from .jit_simulator import QuantumSimulator as JitQuantumSimulator
import time
import sys

type QuantumSimulator = Union[CpuQuantumSimulator, JitQuantumSimulator]

dtype_qstate = np.dtype([("hi", np.uint64), ("lo", np.uint64), ("vec", np.complex128)])


def build_header_from_circuit(circuit: QuantumCircuit) -> dict:
    return {
        "creg_sizes": [[cr.name, cr.size] for cr in circuit.cregs],
        "global_phase": float(circuit.global_phase),
        "memory_slots": sum(cr.size for cr in circuit.cregs),
        "n_qubits": circuit.num_qubits,
        "name": circuit.name,
        "qreg_sizes": [[qr.name, qr.size] for qr in circuit.qregs],
        "metadata": circuit.metadata or {},
    }


def is_zero(val: np.complex128) -> bool:
    return np.abs(val) < 1e-8


def is_one(val: np.complex128) -> bool:
    return np.abs(val - (1.0 + 0j)) < 1e-8


def _x_gate(
    state: NDArray[Any],
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: Tuple[np.uint64, np.uint64],
    ctrl_state: Tuple[np.uint64, np.uint64],
) -> NDArray[Any]:
    bit_hi = np.uint64((1 << target_regs[0]) >> 64)
    bit_lo = np.uint64((1 << target_regs[0]) & 0xFFFFFFFFFFFFFFFF)
    if ctrl_mask[0] or ctrl_mask[1]:
        cond = ((state["hi"] & ctrl_mask[0]) == ctrl_state[0]) & (
            (state["lo"] & ctrl_mask[1]) == ctrl_state[1]
        )
        state["hi"][cond] ^= bit_hi
        state["lo"][cond] ^= bit_lo
    else:
        state["hi"] ^= bit_hi
        state["lo"] ^= bit_lo
    return state


def _diagonal_gate(
    state: NDArray[Any],
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: Tuple[np.uint64, np.uint64],
    ctrl_state: Tuple[np.uint64, np.uint64],
) -> NDArray[Any]:
    bit_hi = np.uint64((1 << target_regs[0]) >> 64)
    bit_lo = np.uint64((1 << target_regs[0]) & 0xFFFFFFFFFFFFFFFF)
    cond0 = ((state["hi"] & bit_hi) == 0) & ((state["lo"] & bit_lo) == 0)
    cond1 = ~cond0

    if ctrl_mask[0] or ctrl_mask[1]:
        cond = ((state["hi"] & ctrl_mask[0]) == ctrl_state[0]) & (
            (state["lo"] & ctrl_mask[1]) == ctrl_state[1]
        )
        cond0 &= cond
        cond1 &= cond
    if not is_one(matrix[0][0]):
        state["vec"][cond0] *= matrix[0][0]
    if not is_one(matrix[1][1]):
        state["vec"][cond1] *= matrix[1][1]
    return state


def _anti_diagonal_gate(
    state: NDArray[Any],
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: Tuple[np.uint64, np.uint64],
    ctrl_state: Tuple[np.uint64, np.uint64],
) -> NDArray[Any]:
    bit_hi = np.uint64((1 << target_regs[0]) >> 64)
    bit_lo = np.uint64((1 << target_regs[0]) & 0xFFFFFFFFFFFFFFFF)
    cond0 = ((state["hi"] & bit_hi) == 0) & ((state["lo"] & bit_lo) == 0)
    cond1 = ~cond0

    if ctrl_mask[0] or ctrl_mask[1]:
        cond = ((state["hi"] & ctrl_mask[0]) == ctrl_state[0]) & (
            (state["lo"] & ctrl_mask[1]) == ctrl_state[1]
        )
        cond0 &= cond
        cond1 &= cond
        state["hi"][cond] ^= bit_hi
        state["lo"][cond] ^= bit_lo
    else:
        state["hi"] ^= bit_hi
        state["lo"] ^= bit_lo

    if not is_one(matrix[1][0]):
        state["vec"][cond0] *= matrix[1][0]
    if not is_one(matrix[0][1]):
        state["vec"][cond1] *= matrix[0][1]
    return state


def _swap_gate(
    state: NDArray[Any],
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: Tuple[np.uint64, np.uint64],
    ctrl_state: Tuple[np.uint64, np.uint64],
) -> NDArray[Any]:
    bit_hi = np.uint64(((1 << target_regs[0]) >> 64) | ((1 << target_regs[1]) >> 64))
    bit_lo = np.uint64(
        ((1 << target_regs[0]) | (1 << target_regs[1])) & 0xFFFFFFFFFFFFFFFF
    )
    cond = ~(
        (((state["hi"] & bit_hi) == 0) & ((state["lo"] & bit_lo) == 0))
        | (((state["hi"] & bit_hi) == bit_hi) & ((state["lo"] & bit_lo) == bit_lo))
    )
    if ctrl_mask[0] or ctrl_mask[1]:
        cond &= ((state["hi"] & ctrl_mask[0]) == ctrl_state[0]) & (
            (state["lo"] & ctrl_mask[1]) == ctrl_state[1]
        )
    state["hi"][cond] ^= bit_hi
    state["lo"][cond] ^= bit_lo
    return state


def _square_gate(
    state: NDArray[Any],
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: Tuple[np.uint64, np.uint64],
    ctrl_state: Tuple[np.uint64, np.uint64],
) -> NDArray[Any]:
    if len(target_regs) > 1:
        raise ValueError("Too many target bits: " + str(len(target_regs)))
    bit_hi = np.uint64((1 << target_regs[0]) >> 64)
    bit_lo = np.uint64((1 << target_regs[0]) & 0xFFFFFFFFFFFFFFFF)

    hi_masked = state["hi"] & ~bit_hi
    lo_masked = state["lo"] & ~bit_lo

    sort_idx = np.lexsort((lo_masked, hi_masked))
    state[:] = state[sort_idx]
    hi_masked[:] = hi_masked[sort_idx]
    lo_masked[:] = lo_masked[sort_idx]

    keep_mask = np.ones(len(state), dtype=bool)
    new_elements = np.empty(len(state), dtype=dtype_qstate)
    next_idx = 0

    i = 0
    while i < len(state):
        if ctrl_mask[0] or ctrl_mask[1]:
            if (hi_masked[i] & ctrl_mask[0]) != ctrl_state[0] or (
                lo_masked[i] & ctrl_mask[1]
            ) != ctrl_state[1]:
                # 制御ビット対象外
                i += 1
                continue
        if (
            i + 1 < len(state)
            and hi_masked[i] == hi_masked[i + 1]
            and lo_masked[i] == lo_masked[i + 1]
        ):
            # ペア
            ix0 = i
            ix1 = i + 1
            if (state["hi"][i] & bit_hi) or (state["lo"][i] & bit_lo):
                ix0, ix1 = ix1, ix0
            vec0 = state["vec"][ix0] * matrix[0][0] + state["vec"][ix1] * matrix[1][0]
            vec1 = state["vec"][ix0] * matrix[0][1] + state["vec"][ix1] * matrix[1][1]
            if is_zero(vec0):
                keep_mask[ix0] = False
            else:
                state["vec"][ix0] = vec0
            if is_zero(vec1):
                keep_mask[ix1] = False
            else:
                state["vec"][ix1] = vec1
            i += 2
        else:
            # シングル
            if (state["hi"][i] & bit_hi) or (state["lo"][i] & bit_lo):
                # 1
                new_elements[next_idx] = (
                    hi_masked[i],
                    lo_masked[i],
                    state["vec"][i] * matrix[0][1],
                )
                state["vec"][i] *= matrix[1][1]
            else:
                # 0
                new_elements[next_idx] = (
                    hi_masked[i] | bit_hi,
                    lo_masked[i] | bit_lo,
                    state["vec"][i] * matrix[1][0],
                )
                state["vec"][i] *= matrix[0][0]
            next_idx += 1
            i += 1

    return np.concatenate([state[keep_mask], new_elements[:next_idx]])


class ResultJob(JobV1):
    """Resultをそのまま返すだけの最小Jobクラス"""

    def __init__(self, backend, result):
        job_id = str(uuid.uuid4())  # 自動で一意なIDを生成
        super().__init__(backend, job_id=job_id)
        self._result = result

    def result(self):
        return self._result

    def submit(self):
        pass

    def status(self):
        # 同期実行なので常にDONE
        return JobStatus.DONE

    def cancel(self):
        # このJobはキャンセルできない
        raise NotImplementedError("This job type cannot be cancelled.")


class SparseStatevectorSimulator(BackendV2):
    # cname, cindex, target_reg
    _measure: List[Tuple[str, int, int]]
    # key=cname, value=measure value
    _result: Dict[str, int]
    # 実行開始時間
    _start_time: float
    # JITモード
    _jit_mode: Union[None, str]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._backend_name = "sparse_simulator"
        self._options = Options(shots=1024)
        self._num_qubits = 90
        self._jit_mode = None

    def _parse_gate(
        self,
        simulator: QuantumSimulator,
        instr: Any,
        target_regs: List[int],
        ctrl_regs: List[int],
        neg_ctrl_regs: List[int],
    ):
        if isinstance(instr, (Measure, Barrier, Reset)):
            # skip
            pass
        elif isinstance(instr, ControlledGate):
            new_ctrl_regs = ctrl_regs[:]
            new_neg_ctrl_regs = neg_ctrl_regs[:]
            for i in range(instr.num_ctrl_qubits):
                if instr.ctrl_state & (1 << i):
                    new_ctrl_regs.append(target_regs[i])
                else:
                    new_neg_ctrl_regs.append(target_regs[i])
            self._parse_gate(
                simulator,
                instr.base_gate,
                target_regs[instr.num_ctrl_qubits :],
                new_ctrl_regs,
                new_neg_ctrl_regs,
            )
        elif hasattr(instr, "definition") and instr.definition is not None:
            try:
                mat = instr.to_matrix()
                simulator.execute(
                    instr.name, mat, target_regs, ctrl_regs, neg_ctrl_regs
                )
            except:
                self._parse_circuit(
                    simulator, instr.definition, target_regs, ctrl_regs, neg_ctrl_regs
                )
        elif hasattr(instr, "to_matrix") and instr.to_matrix is not None:
            simulator.execute(
                instr.name, instr.to_matrix(), target_regs, ctrl_regs, neg_ctrl_regs
            )
        else:
            print("???", instr)

    def _measure_value(self, simulator: QuantumSimulator, m: Tuple[str, int, int]):
        if not m[0] in self._result:
            self._result[m[0]] = 0
        if simulator.measure(m[2]):
            self._result[m[0]] |= 1 << m[1]

    def _measure_count(
        self,
        circuit: QuantumCircuit,
        counts: Dict[str, int],
        memory: Union[List[int], None] = None,
    ):
        keystr = ""
        val = 0
        sht = 0
        for cl in circuit.cregs:
            if cl.name in self._result:
                val |= self._result[cl.name] << sht
            sht += len(cl)
        keystr = hex(val)
        if keystr in counts:
            counts[keystr] += 1
        else:
            counts[keystr] = 1

    def _make_statevector(self, circuit: QuantumCircuit, simulator: QuantumSimulator) -> Generator[Tuple[str, np.complex128]]:
        for st, amp in simulator.get_statevector():
            keystr = ""
            sht = 0
            for qr in circuit.qregs:
                bstr = bin(((st >> sht) & ((1 << len(qr)) - 1)))[2:].zfill(len(qr))
                if keystr:
                    keystr = " " + keystr
                keystr = bstr + keystr
                sht += len(qr)
            yield (keystr, amp)

    def _parse_circuit(
        self,
        simulator: QuantumSimulator,
        circuit: QuantumCircuit,
        target_regs: Union[List[int], None] = None,
        ctrl_regs: Union[List[int], None] = None,
        neg_ctrl_regs: Union[List[int], None] = None,
    ):
        self._start_time = time.perf_counter()
        debug_out = False
        if target_regs is None:
            target_regs = list(range(len(circuit.qubits)))
            debug_out = True
        ix = 0
        size = len(circuit.data)
        for instr, qargs, cargs in circuit.data:
            ix += 1
            if debug_out and (ix % 10000) == 0:
                tm = time.perf_counter() - self._start_time
                hour = int(tm // 3600)
                min = int((tm % 3600) // 60)
                sec = tm % 60
                if hour > 0:
                    tmstr = f"{hour:d}:{min:02d}:{sec:06.3f}"
                elif min > 0:
                    tmstr = f"{min:d}:{sec:06.3f}"
                else:
                    tmstr = f"{sec:.3f}"
                print(f"execute circuit {ix}/{size} {tmstr}")
            if isinstance(instr, Measure):
                for ix in range(len(qargs)):
                    self._measure.append(
                        (
                            cargs[ix]._register.name,
                            cargs[ix]._index,
                            target_regs[circuit.qubits.index(qargs[ix])],
                        )
                    )
                    self._measure_value(simulator, self._measure[-1])
                continue
            elif isinstance(instr, Reset):
                for ix in range(len(qargs)):
                    simulator.reset(target_regs[circuit.qubits.index(qargs[ix])])
                continue
            self._parse_gate(
                simulator,
                instr,
                [target_regs[circuit.qubits.index(qb)] for qb in qargs],
                ctrl_regs or [],
                neg_ctrl_regs or [],
            )

    @property
    def target(self):
        pass
        # return self._target

    @property
    def num_qubits(self):
        return self._num_qubits

    @property
    def max_circuits(self):
        return 1

    @property
    def options(self):
        return self._options

    @classmethod
    def _default_options(cls):
        return Options(shots=1024)

    def set_jit_mode(self, mode: Union[None, str]):
        """JITモードを設定する。mode=None, "single", "parallel" のいずれか"""
        self._jit_mode = mode

    def run(self, circuits, **run_options) -> JobV1:
        """Qiskitが実際に回路を実行する際に呼び出すメソッド"""
        if not isinstance(circuits, list):
            circuits = [circuits]

        self._measure = []
        self._result = {}

        if self._jit_mode == "single":
            simulator = JitQuantumSimulator(circuits[0].num_qubits)
        elif self._jit_mode == "parallel":
            simulator = JitQuantumSimulator(circuits[0].num_qubits, True)
        else:
            simulator = CpuQuantumSimulator(circuits[0].num_qubits)
        self._parse_circuit(simulator, circuits[0])  # type: ignore

        # simulator.dump()
        counts = {}
        mem_list = None
        if run_options.get("memory"):
            mem_list = []
        self._measure_count(circuits[0], counts, mem_list)
        shots = run_options.get("shots", 1)
        for i in range(1, shots):
            simulator.reset_measure()
            self._result = {}
            for mes in self._measure:
                self._measure_value(simulator, mes)  # type: ignore
            self._measure_count(circuits[0], counts, mem_list)

        mydt = ExperimentResultData(counts=counts, memory=mem_list)
        myres = ExperimentResult(
            shots, True, mydt, header=build_header_from_circuit(circuits[0])
        )
        result = Result(results=[myres], quasi_dists=[{"01": 0.5}], state_vector=self._make_statevector(circuits[0], simulator))

        return ResultJob(self, result)


def get_backend(
    circuit: QuantumCircuit,
) -> Tuple[SparseStatevectorSimulator, QuantumCircuit]:
    jit_mode = None
    for arg in sys.argv:
        if arg == "-single":
            jit_mode = "single"
        elif arg == "-parallel":
            jit_mode = "parallel"
    simulator = SparseStatevectorSimulator()
    print(f"JIT mode: {jit_mode}")
    simulator.set_jit_mode(jit_mode)
    return (simulator, circuit)
