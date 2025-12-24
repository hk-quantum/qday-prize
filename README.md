# QDays Prize

## Email Address

hk.quantum@icloud.com

## Background

I am a software engineer at a Japanese software company. To better understand how quantum computers work, I developed my own quantum circuit simulator and tooling.

For RSA cryptanalysis (factoring composite numbers made from two primes), I have simulated quantum circuits that succeeded for these bit sizes using my custom tools:

- Broke 10 bits with Shor's algorithm
- Broke 40 bits with Grover's algorithm


## What I Broke This Time

The largest bit size for which I was able to run a quantum circuit on the simulator is 12 bits:

```
--- Bit size 12 ---
Bit size: 12
Prime p: 2089
Curve order (#E): 2143
Subgroup order n: 2143
Cofactor h: 1
Generator point G: (1417, 50)
Private key d: 1384
Public key Q: (1043, 1795)
```

The largest bit size executed on an IBM quantum computer (ibm_fez) in this work was 5 bits:

```
"bit_length": 5,
"prime": 23,
"a": 1,
"b": 11,
"generator_point": [7, 4],
"public_key": [9, 6]
```


## Execution Environment

Quantum circuits are implemented using Qiskit and executed with my high-performance simulator `SparseStatevectorSimulator` as the backend.

On IBM's quantum device `ibm_fez`, I could run a 5-bit ECC curve, but due to circuit depth and device noise the results were essentially random and did not reach the expected accuracy. Measured execution times for 100 shots were approximately:

- 3 bits: ~3 seconds
- 4 bits: ~8 seconds
- 5 bits: ~30 seconds


## Execution Procedure

### Initial Setup

Create a Python virtual environment and install dependencies:

```
git clone https://github.com/hk-quantum/qday-prize.git
cd qday-prize
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

To run on IBM Quantum hardware, set your IBM Quantum Platform API key and CRN in a `.env` file:

```
API_KEY=<Your API Key>
IBM_CRN=<Your CRN>
```


### How to Run

Specify the bit size to target and the algorithm variant (`compact` or `wide`) as command-line arguments. Optionally choose single-threaded (`-single`) or multi-threaded (`-parallel`) JIT execution.

```
python src/main.py <num_bits> compact|wide [-single|-parallel]
```

The program reads ECC parameters for the chosen bit size from `data/curves.json`, constructs the quantum circuit, runs the selected backend, and writes execution logs. At the end it reports the number of successful private-key recoveries out of 100 measurement shots.

To run on IBM hardware use:

```
python src/ibm_main.py <num_bits> compact|wide
```


### Execution Results

From the 100-shot measurements the recovered private key and number of successful recoveries are printed in the last line of the log.

Example runs on IBM hardware (`ibm_fez`):

| log file | type    | output result |
|---|---|---|
| `logs/3_ibm_compact.out` | compact | `Success: d=3 count=33` |
| `logs/3_ibm_wide.out`    | wide    | `Success: d=3 count=34` |
| `logs/4_ibm_compact.out` | compact | `Success: d=6 count=27` |
| `logs/4_ibm_wide.out`    | wide    | `Success: d=6 count=26` |
| `logs/5_ibm_compact.out` | compact | `Success: d=6 count=24` |


The following results were obtained using my simulator `SparseStatevectorSimulator` (see [Version 4](https://github.com/hk-quantum/qday-prize/tree/v4.0.0)):

| log file | type    | output result |
|---|---|---|
| `logs/3_compact.out`  | compact | `Success: d=3 count=68` |
| `logs/3_wide.out`     | wide    | `Success: d=3 count=75` |
| `logs/4_compact.out`  | compact | `Success: d=6 count=79` |
| `logs/4_wide.out`     | wide    | `Success: d=6 count=75` |
| `logs/5_compact.out`  | compact | `Success: d=6 count=80` |
| `logs/5_wide.out`     | wide    | `Success: d=6 count=75` |
| `logs/6_compact.out`  | compact | `Success: d=18 count=82` |
| `logs/6_wide.out`     | wide    | `Success: d=18 count=76` |
| `logs/7_compact.out`  | compact | `Success: d=56 count=85` |
| `logs/7_wide.out`     | wide    | `Success: d=56 count=85` |
| `logs/8_compact.out`  | compact | `Success: d=103 count=87` |
| `logs/8_wide.out`     | wide    | `Success: d=103 count=88` |
| `logs/9_compact.out`  | compact | `Success: d=135 count=87` |
| `logs/9_wide.out`     | wide    | `Success: d=135 count=85` |
| `logs/10_compact.out` | compact | `Success: d=165 count=88` |
| `logs/10_wide.out`    | wide    | `Success: d=165 count=95` |
| `logs/11_compact.out` | compact | `Success: d=756 count=92` |
| `logs/12_compact.out` | compact | `Success: d=1384 count=88` |

For larger bit sizes, only the `compact` variant is shown where the simulator reached its performance limits.
