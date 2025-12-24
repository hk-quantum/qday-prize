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

from dotenv import load_dotenv
import os
from qiskit import transpile, QuantumCircuit
from qiskit_ibm_runtime import QiskitRuntimeService, Options
from qiskit_ibm_runtime import SamplerV2 as Sampler
from typing import Tuple
import time
import main

load_dotenv()

API_KEY = os.environ["API_KEY"]
IBM_CRN = os.environ["IBM_CRN"]


def get_backend(circuit: QuantumCircuit) -> Tuple[Sampler, QuantumCircuit]:
    QiskitRuntimeService.save_account(
        token=API_KEY,
        instance=IBM_CRN,
        plans_preference=["open"],
        region="us-east",
        set_as_default=True,
        overwrite=True,
    )
    service = QiskitRuntimeService()
    backend = None
    for bk in service.backends(min_num_qubits=circuit.num_qubits):
        if backend is None:
            backend = bk
        elif bk.status().pending_jobs < backend.status().pending_jobs:
            backend = bk
    if backend is None:
        raise RuntimeError("Not Found Backend")

    start_time = time.perf_counter()
    print("Start Transpile", backend.name, backend.status().pending_jobs)
    t_qc = transpile(circuit, backend)
    tm = time.perf_counter() - start_time
    print(f"End Transpile {tm:.3f}[s] size={t_qc.size()} depth={t_qc.depth()}, ops={t_qc.count_ops()}")
    return Sampler(mode=backend), t_qc

if __name__ == "__main__":
    main.main(get_backend)