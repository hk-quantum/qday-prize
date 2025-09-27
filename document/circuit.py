from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.result import Result
from qiskit.result.models import ExperimentResultData, ExperimentResult
from qiskit.circuit import Gate


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


def print_overview():
    reg_a = QuantumRegister(1, name="qa")
    reg_b = QuantumRegister(1, name="qb")
    reg_x = QuantumRegister(2, name="qx")
    reg_ext = QuantumRegister(4, name="anc")
    cl_a = ClassicalRegister(len(reg_a), name="a")
    cl_b = ClassicalRegister(len(reg_b), name="b")
    cl_xy = ClassicalRegister(len(reg_x), name="xy")
    qc = QuantumCircuit(reg_a, reg_b, reg_x, reg_ext, cl_a, cl_b, cl_xy)
    qc.h(reg_a)
    qc.h(reg_b)
    qc.append(
        QuantumCircuit(len(reg_x) + len(reg_ext), name="oracle: aG+bQ")
        .to_gate()
        .control(len(reg_a) + len(reg_b)),
        [*reg_a, *reg_b, *reg_x, *reg_ext],
    )
    qc.measure(reg_x, cl_xy)
    qc.barrier()
    qc.append(QuantumCircuit(len(reg_a), name="QFT†").to_gate(), reg_a)
    qc.append(QuantumCircuit(len(reg_b), name="QFT†").to_gate(), reg_b)
    qc.measure(reg_a, cl_a)
    qc.measure(reg_b, cl_b)
    print(qc)


def print_addr_mod():
    reg_x = QuantumRegister(3, name="qx")
    reg_z = QuantumRegister(3, name="qz")
    reg_ext = QuantumRegister(1, name="anc")
    qc = QuantumCircuit(reg_x, reg_z, reg_ext)
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="z+=x").to_gate().control(len(reg_x)),
        [*reg_x, *reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="z-=n"),
        [*reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z), name="z+=n (if anc)").to_gate().control(1),
        [*reg_ext, *reg_z],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="z-=x").to_gate().control(len(reg_x)),
        [*reg_x, *reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z), name="z+=x").to_gate().control(len(reg_x)),
        [*reg_x, *reg_z],
    )
    qc.x(reg_ext)
    print(qc)


def print_addr_mod_const():
    reg_z = QuantumRegister(3, name="qz")
    reg_ext = QuantumRegister(1, name="anc")
    qc = QuantumCircuit(reg_z, reg_ext)
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="z+=a").to_gate(),
        [*reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="z-=n"),
        [*reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z), name="z+=n (if anc)").to_gate().control(1),
        [*reg_ext, *reg_z],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="z-=a").to_gate(),
        [*reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z), name="z+=a").to_gate(),
        [*reg_z],
    )
    qc.x(reg_ext)
    print(qc)


def print_multi_mod():
    reg_x = QuantumRegister(3, "qx")
    reg_y = QuantumRegister(3, "qy")
    reg_z = QuantumRegister(3, "qz")
    reg_ext = QuantumRegister(1, "anc")
    qc = QuantumCircuit(reg_x, reg_y, reg_z, reg_ext)
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="addr_mod(2^0%n,n)").to_gate().control(2),
        [reg_x[0], reg_y[0], *reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="addr_mod(2^1%n,n)").to_gate().control(2),
        [reg_x[1], reg_y[0], *reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="addr_mod(2^2%n,n)").to_gate().control(2),
        [reg_x[2], reg_y[0], *reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="addr_mod(2^1%n,n)").to_gate().control(2),
        [reg_x[0], reg_y[1], *reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="addr_mod(2^2%n,n)").to_gate().control(2),
        [reg_x[1], reg_y[1], *reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="addr_mod(2^3%n,n)").to_gate().control(2),
        [reg_x[2], reg_y[1], *reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="addr_mod(2^2%n,n)").to_gate().control(2),
        [reg_x[0], reg_y[2], *reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="addr_mod(2^3%n,n)").to_gate().control(2),
        [reg_x[1], reg_y[2], *reg_z, *reg_ext],
    )
    qc.append(
        QuantumCircuit(len(reg_z) + 1, name="addr_mod(2^4%n,n)").to_gate().control(2),
        [reg_x[2], reg_y[2], *reg_z, *reg_ext],
    )
    print(qc)


def print_xy_add_P():
    reg_x = QuantumRegister(1, name="qx")
    reg_y = QuantumRegister(1, name="qy")
    reg_ox = QuantumRegister(1, name="ox")
    reg_oy = QuantumRegister(1, name="oy")
    reg_anc = QuantumRegister(3, name="anc")
    qc = QuantumCircuit(reg_x, reg_y, reg_ox, reg_oy, reg_anc)
    qc.append(
        QuantumCircuit(len(reg_ox) + len(reg_oy) + len(reg_anc), name="(qx,qy)+P")
        .to_gate()
        .control(len(reg_x) + len(reg_y)),
        [*reg_x, *reg_y, *reg_ox, *reg_oy, *reg_anc],
    )
    qc.swap(reg_x, reg_ox)
    qc.swap(reg_y, reg_oy)
    qc.append(
        QuantumCircuit(len(reg_ox) + len(reg_oy) + len(reg_anc), name="inv (qx,qy)-P")
        .to_gate()
        .control(len(reg_x) + len(reg_y)),
        [*reg_x, *reg_y, *reg_ox, *reg_oy, *reg_anc],
    )
    print(qc)


def print_xy_add_P_calc_lambda():
    reg_x = QuantumRegister(1, name="qx")
    reg_y = QuantumRegister(1, name="qy")
    reg_ox = QuantumRegister(1, name="ox")
    reg_oy = QuantumRegister(1, name="oy")
    reg_dx = QuantumRegister(1, name="dx")
    reg_dx_dg = QuantumRegister(1, name="dx^-1")
    reg_dy = QuantumRegister(1, name="dy")
    reg_lm = QuantumRegister(1, name="lambda")
    reg_anc = QuantumRegister(3, name="anc")
    qc = QuantumCircuit(
        reg_x, reg_y, reg_ox, reg_oy, reg_dx, reg_dx_dg, reg_dy, reg_lm, reg_anc
    )
    qc.append(
        QuantumCircuit(len(reg_ox) + 1, name="qx-Px -> dx")
        .to_gate()
        .control(len(reg_x)),
        [*reg_x, *reg_dx, reg_dy[0]],
    )
    qc.append(Gate("dx", len(reg_dx), params=[]), reg_dx)
    qc.append(
        QuantumCircuit(
            len(reg_dx_dg) + len(reg_dy) + len(reg_anc) + len(reg_lm),
            name="x^(p-2)",
        )
        .to_gate()
        .control(len(reg_dx)),
        [*reg_dx, *reg_dx_dg, *reg_dy, *reg_lm, *reg_anc],
    )
    qc.append(Gate("dx^-1", len(reg_dx), params=[]), reg_dx_dg)
    qc.append(
        QuantumCircuit(len(reg_ox) + 1, name="qy-Py -> dy")
        .to_gate()
        .control(len(reg_x)),
        [*reg_y, *reg_dy, reg_lm[0]],
    )
    qc.append(Gate("dy", len(reg_dy), params=[]), reg_dy)
    qc.append(
        QuantumCircuit(1, name="(qx,qy)==P").to_gate().control(len(reg_x) + len(reg_y)),
        [*reg_x, *reg_y, reg_anc[1]],
    )
    qc.append(
        QuantumCircuit(len(reg_lm) + 1, name="dx^-1 * dy mod n")
        .to_gate()
        .control(
            len(reg_dx_dg) + len(reg_dy) + 1,
            ctrl_state=(1 << (len(reg_dx_dg) + len(reg_dy))) - 1,
        ),
        [*reg_dx_dg, *reg_dy, reg_anc[1], *reg_lm, reg_anc[0]],
    )
    qc.append(
        QuantumCircuit(len(reg_lm), name="(3*Px^2+a)/(2*Py)").to_gate().control(1),
        [reg_anc[1], *reg_lm],
    )
    qc.append(Gate("lambda", len(reg_lm), params=[]), reg_lm)
    qc.barrier()

    print(qc)


def print_xy_add_P_with_lambda():
    reg_x = QuantumRegister(1, name="qx")
    reg_y = QuantumRegister(1, name="qy")
    reg_ox = QuantumRegister(1, name="ox")
    reg_oy = QuantumRegister(1, name="oy")
    reg_dx = QuantumRegister(1, name="dx")
    reg_dx_dg = QuantumRegister(1, name="dx^-1")
    reg_dy = QuantumRegister(1, name="dy")
    reg_lm = QuantumRegister(1, name="lambda")
    reg_anc = QuantumRegister(3, name="anc")
    qc = QuantumCircuit(
        reg_x, reg_y, reg_ox, reg_oy, reg_dx, reg_dx_dg, reg_dy, reg_lm, reg_anc
    )
    qc.append(
        QuantumCircuit(len(reg_ox) + 1, name="lambda^2").to_gate().control(len(reg_lm)),
        [*reg_lm, *reg_ox, reg_oy[0]],
    )
    qc.append(
        QuantumCircuit(len(reg_ox) + 1, name="ox-Px").to_gate(), [*reg_ox, reg_oy[0]]
    )
    qc.append(
        QuantumCircuit(len(reg_ox) + 1, name="ox-x").to_gate().control(len(reg_x)),
        [*reg_x, *reg_ox, reg_oy[0]],
    )
    qc.append(Gate("ox", len(reg_oy), params=[]), reg_ox)
    qc.append(
        QuantumCircuit(len(reg_oy) + 1, name="lmbda * ox")
        .to_gate()
        .control(len(reg_ox) + len(reg_lm)),
        [*reg_ox, *reg_lm, *reg_oy, reg_anc[0]],
    )
    qc.append(
        QuantumCircuit(2, name="lmbda * P").to_gate().control(len(reg_lm)),
        [*reg_lm, *reg_anc[:2]],
    )
    print(qc)


# print_overview()
# print_addr_mod()
# print_addr_mod_const()
# print_multi_mod()
#print_xy_add_P_with_lambda()
print_xy_add_P()