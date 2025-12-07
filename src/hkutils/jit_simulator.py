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
from typing import List, Tuple, Any, Union, Dict, Callable, Generator, Literal
from numpy.typing import NDArray
import time
import math
from collections import defaultdict
from numba import njit, prange, uint32, void


def is_zero(val: np.complex128) -> bool:
    return np.abs(val) < 1e-8


def is_one(val: np.complex128) -> bool:
    return np.abs(val - (1.0 + 0j)) < 1e-8


def get_ctrl_mask(ctrl_regs: List[int], neg_ctrl_regs: List[int]) -> NDArray[np.uint32]:
    """
    ctrl_regs: 正制御ビット位置リスト
    neg_ctrl_regs: 負制御ビット位置リスト
    -> (index, mask, state) のリスト（indexごとに合成）
    """

    mask_dict = defaultdict(lambda: [0, 0])  # index: [mask, state]
    for reg in ctrl_regs:
        idx, bit = divmod(reg, 32)
        mask_dict[idx][0] |= 1 << bit
        mask_dict[idx][1] |= 1 << bit  # 正制御はstate=1
    for reg in neg_ctrl_regs:
        idx, bit = divmod(reg, 32)
        mask_dict[idx][0] |= 1 << bit
        # 負制御はstate=0（stateの該当bitは0のまま）
    result = np.array(
        [(idx, mask, state) for idx, (mask, state) in mask_dict.items() if mask != 0],
        dtype=np.uint32,
    )
    return result

@njit(parallel=True)
def _x_gate_jit_parallel(
    qstate: NDArray[np.uint32],  # dtype_qstate の配列
    matrix: NDArray[np.complex128],
    iy: np.uint32,
    bit: np.uint32,
    ctrl_mask: NDArray[np.uint32],
) -> np.ndarray:
    """
    数千万行のqstateに対して条件付きXゲートを高速に適用する。
    qstate["state"] は shape=(N, 8) の np.uint32 配列を想定 (256bit状態)
    """
    for i in prange(qstate.shape[0]):
        # --- コントロール条件判定 ---
        cond = True
        for cm in ctrl_mask:
            idx = cm[0]
            mask = cm[1]
            cstate = cm[2]
            if (qstate[i, idx] & mask) != cstate:
                cond = False
                break
        if not cond:
            continue

        # --- 条件成立行のみXOR ---
        qstate[i, iy] ^= bit

    return qstate

@njit(parallel=False)
def _x_gate_jit(
    qstate: NDArray[np.uint32],  # dtype_qstate の配列
    matrix: NDArray[np.complex128],
    iy: np.uint32,
    bit: np.uint32,
    ctrl_mask: NDArray[np.uint32],
) -> np.ndarray:
    """
    数千万行のqstateに対して条件付きXゲートを高速に適用する。
    qstate["state"] は shape=(N, 8) の np.uint32 配列を想定 (256bit状態)
    """
    for i in prange(qstate.shape[0]):
        # --- コントロール条件判定 ---
        cond = True
        for cm in ctrl_mask:
            idx = cm[0]
            mask = cm[1]
            cstate = cm[2]
            if (qstate[i, idx] & mask) != cstate:
                cond = False
                break
        if not cond:
            continue

        # --- 条件成立行のみXOR ---
        qstate[i, iy] ^= bit

    return qstate

def _x_gate(
    qstate: np.ndarray,  # dtype_qstate の配列
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: NDArray[np.uint32],
) -> np.ndarray:
    state = qstate["state"]

    # --- ターゲットビット計算 ---
    iy, ix = divmod(target_regs[0], 32)
    bit = np.uint32(1 << ix)

    # --- コントロールなしパス: 最速処理 ---
    if ctrl_mask.shape[0] == 0:
        state[:, iy] ^= bit
        return qstate
    
    _x_gate_jit(state, matrix, np.uint32(iy), np.uint32(bit), ctrl_mask)
    return qstate

def _x_gate_parallel(
    qstate: np.ndarray,  # dtype_qstate の配列
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: NDArray[np.uint32],
) -> np.ndarray:
    state = qstate["state"]

    # --- ターゲットビット計算 ---
    iy, ix = divmod(target_regs[0], 32)
    bit = np.uint32(1 << ix)

    # --- コントロールなしパス: 最速処理 ---
    if ctrl_mask.shape[0] == 0:
        state[:, iy] ^= bit
        return qstate
    
    _x_gate_jit_parallel(state, matrix, np.uint32(iy), np.uint32(bit), ctrl_mask)
    return qstate


def _diagonal_gate(
    qstate: np.ndarray,  # dtype_qstate の配列
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: NDArray[np.uint32],
) -> NDArray[Any]:
    """
    target_regs のビットが 0/1 に応じて vec を diagonal matrix の対角成分で掛ける
    ctrl_mask が指定されていれば、その条件を満たす行のみ適用
    """
    iy, ix = divmod(target_regs[0], 32)
    bit = np.uint32(1 << ix)

    apply0 = (qstate["state"][:, iy] & bit) == 0
    apply1 = ~apply0

    # ctrl_mask がある場合は mask 計算
    if ctrl_mask.shape[0] > 0:
        mask = np.ones((len(qstate),), dtype=bool)
        for cm in ctrl_mask:
            idx, mask_val, cstate = cm
            mask &= (qstate["state"][:, idx] & mask_val) == cstate
        apply0 &= mask
        apply1 &= mask

    if not is_one(matrix[0, 0]):
        qstate["vec"][apply0] *= matrix[0, 0]
    if not is_one(matrix[1, 1]):
        qstate["vec"][apply1] *= matrix[1, 1]

    return qstate


def _anti_diagonal_gate(
    qstate: np.ndarray,  # dtype_qstate の配列
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: NDArray[np.uint32],
) -> NDArray[Any]:
    iy, ix = divmod(target_regs[0], 32)
    bit = np.uint32(1 << ix)

    apply0 = (qstate["state"][:, iy] & bit) == 0
    apply1 = ~apply0

    # ctrl_mask がある場合は mask 計算
    if ctrl_mask.shape[0] > 0:
        mask = np.ones((len(qstate),), dtype=bool)
        for cm in ctrl_mask:
            idx, mask_val, cstate = cm
            mask &= (qstate["state"][:, idx] & mask_val) == cstate
        apply0 &= mask
        apply1 &= mask
        qstate["state"][mask, iy] ^= bit
    else:
        qstate["state"][:, iy] ^= bit

    if not is_one(matrix[1, 0]):
        qstate["vec"][apply0] *= matrix[1, 0]
    if not is_one(matrix[0, 1]):
        qstate["vec"][apply1] *= matrix[0, 1]

    return qstate


def _swap_gate(
    qstate: np.ndarray,  # dtype_qstate の配列
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: NDArray[np.uint32],
) -> NDArray[Any]:
    """
    target_regs[0] と target_regs[1] のいずれか片方だけが1の場合、
    両方のビットを反転させる。
    ctrl_mask があれば条件付きで適用
    """
    N = qstate.shape[0]

    # ターゲットビットの位置
    iy0, ix0 = divmod(target_regs[0], 32)
    iy1, ix1 = divmod(target_regs[1], 32)
    bit0 = np.uint32(1 << ix0)
    bit1 = np.uint32(1 << ix1)

    # まず XOR 条件を計算
    b0 = (qstate["state"][:, iy0] & bit0) != 0
    b1 = (qstate["state"][:, iy1] & bit1) != 0
    mask = b0 != b1  # どちらか片方だけ 1 の行

    # ctrl_mask 条件がある場合はマージ
    if ctrl_mask.shape[0] > 0:
        for cm in ctrl_mask:
            idx, mask_val, cstate = cm
            mask &= (qstate["state"][:, idx] & mask_val) == cstate

    # 条件付きで両方のビットを反転
    qstate["state"][mask, iy0] ^= bit0
    qstate["state"][mask, iy1] ^= bit1

    return qstate


def _square_gate(
    qstate: np.ndarray,  # dtype_qstate の配列
    matrix: NDArray[np.complex128],
    target_regs: List[int],
    ctrl_mask: NDArray[np.uint32],
) -> np.ndarray:
    if len(target_regs) > 1:
        raise ValueError("Too many target bits: " + str(len(target_regs)))

    # ターゲットビット位置
    iy, ix = divmod(target_regs[0], 32)
    bit = np.uint32(1 << ix)

    # state_masked = 条件付き XOR 用のコピー（ソート用）
    state_masked = qstate["state"].copy()
    state_masked[:, iy] &= ~bit  # 対象ビットを 0 にしてソート

    # ソートしてペアを作りやすくする
    keys = [state_masked[:, i] for i in reversed(range(state_masked.shape[1]))]
    sort_idx = np.lexsort(keys)
    state = qstate[sort_idx]
    state_masked = state_masked[sort_idx]

    keep_mask = np.ones(len(state), dtype=bool)
    new_elements = np.empty(len(state), dtype=state.dtype)
    next_idx = 0

    i = 0
    while i < len(state):
        # ctrl_mask 条件を満たさなければスキップ
        if ctrl_mask.shape[0] > 0:
            ng = False
            for cm in ctrl_mask:
                idx, mask_val, cstate = cm
                if (state["state"][i, idx] & mask_val) != cstate:
                    ng = True
                    break
            if ng:
                i += 1
                continue

        # ペア判定
        if i + 1 < len(state) and np.array_equal(state_masked[i], state_masked[i + 1]):
            ix0, ix1 = i, i + 1
            if state["state"][i, iy] & bit:
                ix0, ix1 = ix1, ix0

            vec0 = state["vec"][ix0] * matrix[0, 0] + state["vec"][ix1] * matrix[1, 0]
            vec1 = state["vec"][ix0] * matrix[0, 1] + state["vec"][ix1] * matrix[1, 1]

            if is_zero(vec0):  # 以前の is_zero 相当
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
            if state["state"][i, iy] & bit:
                # 1 の場合
                new_elements[next_idx]["state"] = state["state"][i]
                new_elements[next_idx]["state"][iy] &= ~bit
                new_elements[next_idx]["vec"] = state["vec"][i] * matrix[0, 1]
                state["vec"][i] *= matrix[1, 1]
            else:
                # 0 の場合
                new_elements[next_idx]["state"] = state["state"][i].copy()
                new_elements[next_idx]["state"][iy] |= bit
                new_elements[next_idx]["vec"] = state["vec"][i] * matrix[1, 0]
                state["vec"][i] *= matrix[0, 0]

            next_idx += 1
            i += 1

    return np.ascontiguousarray(np.concatenate([state[keep_mask], new_elements[:next_idx]]))

class QuantumSimulator:
    _state: np.ndarray
    _def_gete: Dict[
        str,
        Callable[
            [
                np.ndarray,
                NDArray[np.complex128],
                List[int],
                NDArray[np.uint32],
            ],
            np.ndarray,
        ],
    ]
    _measure_mask: int = 0
    _measure_state: int = 0
    _max_qubits = 0

    def __init__(self, bit_num, parallel: bool = False):
        # self._max_qubits = bit_num
        num = math.ceil(bit_num / 32)
        self._state = np.array(
            [(np.zeros([num], dtype=np.uint32), 1 + 0j)],
            dtype=np.dtype([("state", np.uint32, num), ("vec", np.complex128)]),
        )
        if parallel:
            self._def_gete = {"x": _x_gate_parallel, "swap": _swap_gate}
        else:
            self._def_gete = {"x": _x_gate, "swap": _swap_gate}

    def execute(
        self,
        name: str,
        matrix: NDArray[np.complex128],
        target_regs: List[int],
        ctrl_regs: List[int],
        neg_ctrl_regs: List[int],
    ):
        for qb in target_regs:
            self._max_qubits = max(self._max_qubits, qb + 1)
        self.apply_measure()
        ctrl_mask = get_ctrl_mask(ctrl_regs, neg_ctrl_regs)
        try:
            if name in self._def_gete:
                self._state = self._def_gete[name](
                    self._state, matrix, target_regs, ctrl_mask
                )
            elif len(target_regs) == 1:
                proc = _square_gate
                if is_zero(matrix[1][0]) and is_zero(matrix[0][1]):
                    proc = _diagonal_gate
                elif is_zero(matrix[0][0]) and is_zero(matrix[1][1]):
                    proc = _anti_diagonal_gate
                self._state = proc(
                    self._state, matrix, target_regs, ctrl_mask
                )
            else:
                raise ValueError("Not Support multi bit gate")
        except:
            import traceback

            traceback.print_exc()

    def measure(self, target: int) -> int:
        bit_flag = 1 << target
        iy, ix = divmod(target, 32)
        bit = np.uint32(1 << ix)

        if (self._measure_mask & bit_flag) == 0:
            # ターゲットビットの条件
            cond0 = (self._state["state"][:, iy] & bit) == 0
            cond1 = ~cond0

            # measure_mask がある場合はマージ
            if self._measure_mask:
                mask_cond = np.ones((len(self._state),), dtype=bool)
                tmp_mask = self._measure_mask
                tmp_state = self._measure_state
                fx = 0
                while tmp_mask > 0:
                    if tmp_mask & 0xFFFFFFFF:
                        mask_cond &= (
                            self._state["state"][:, fx]
                            & np.uint32(tmp_mask & 0xFFFFFFFF)
                        ) == (tmp_state & 0xFFFFFFFF)
                    tmp_mask >>= 32
                    tmp_state >>= 32
                    fx += 1
                cond0 &= mask_cond
                cond1 &= mask_cond

            # 確率計算
            sum0 = np.sum(np.abs(self._state["vec"][cond0]) ** 2)
            sum1 = np.sum(np.abs(self._state["vec"][cond1]) ** 2)

            # 測定結果をランダムに決定
            total = float(sum0 + sum1)
            r = np.random.uniform(0, total)
            if r >= sum0:
                self._measure_state |= bit_flag

            # measure_mask を更新
            self._measure_mask |= bit_flag

        return 1 if self._measure_state & bit_flag else 0

    def dump(self):
        sort_idx = np.argsort(-np.abs(self._state["vec"]))
        self._state = self._state[sort_idx]
        for i in range(len(self._state)):
            st = ""
            for j in range(0, self._max_qubits, 32):
                bit_len = min(32, self._max_qubits - j)
                st = format(self._state["state"][i, j // 32], f"0{bit_len}b") + st
            print(
                st,
                self._state["vec"][i],
            )
        print("Qubit Num", self._max_qubits, "StateCount", len(self._state))

    def get_statevector(self, sort_by: Union[Literal["state", "vec"], None]=None) -> Generator[Tuple[int, np.complex128]]:
        if sort_by == "state":
            L = self._state["state"].shape[0]
            # lexsort のキーは「逆順に」渡す
            keys = [self._state["state"][:, i] for i in range(L-1, -1, -1)]
            sort_idx = np.lexsort(keys)
            self._state = self._state[sort_idx]
        elif sort_by == "vec":
            sort_idx = np.argsort(-np.abs(self._state["vec"]))
            self._state = self._state[sort_idx]
        for i in range(len(self._state)):
            st = 0
            shft = 0
            for v in self._state["state"][i]:
                st |= v << shft
                shft += 32
            yield (st, self._state["vec"][i])

    def reset(self, target: int):
        val = self.measure(target)
        if val:
            self.execute(
                "x",
                np.array([[0j, 1 + 0j], [1 + 0j, 0j]], dtype=np.complex128),
                [target],
                [],
                [],
            )

    def reset_measure(self):
        self._measure_mask = 0
        self._measure_state = 0

    def apply_measure(self):
        if self._measure_mask == 0:
            return

        cond = np.ones((len(self._state),), dtype=bool)
        tmp_mask = self._measure_mask
        tmp_state = self._measure_state
        fx = 0
        while tmp_mask > 0:
            if tmp_mask & 0xFFFFFFFF:
                cond &= (
                    self._state["state"][:, fx] & np.uint32(tmp_mask & 0xFFFFFFFF)
                ) == (tmp_state & 0xFFFFFFFF)
            tmp_mask >>= 32
            tmp_state >>= 32
            fx += 1

        # 条件を満たす行だけ残す
        self._state = self._state[cond]

        # vec を正規化
        total = np.sum(np.abs(self._state["vec"] ** 2))
        if total > 0:
            self._state["vec"] /= np.sqrt(total)

        # measure_mask/state をリセット
        self._measure_mask = 0
        self._measure_state = 0
