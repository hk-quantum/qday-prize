# QDays Prize

## Email Address

hk.quantum@icloud.com

## Background

I am a software developer at a Japanese software development company. To understand how quantum computers work, I created my own quantum computer simulator and related tools.

For RSA cryptanalysis (factoring a composite number made from two primes), I have succeeded in simulating quantum circuits that break the following bit sizes using my custom simulator:

- Broke 10 bits with Shor's algorithm
- Broke 40 bits with Grover's algorithm


## What I Broke This Time

The maximum bit size for which I was able to execute a quantum circuit on the simulator was 12 bits.

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

The maximum bit size executed on a real IBM quantum computer (ibm_torino) was 4 bits.

```
--- Bit size 4 ---
Bit size: 4
Prime p: 13
Curve order (#E): 7
Subgroup order n: 7
Cofactor h: 1
Generator point G: (11, 5)
Private key d: 6
Public key Q: (11, 8)
```


## Execution Environment

The quantum circuits are implemented with Qiskit and executed using my high-performance custom simulator `SparseStatevectorSimulator` as the backend.

On IBM's quantum computer "ibm_torino", I was able to run a 4-bit ECC curve, but due to circuit depth and noise the results were essentially random and did not achieve the expected accuracy. The execution time for 100 shots was about 35 seconds.

## Execution Procedure


### Initial Setup

Set up a Python virtual environment and install the required libraries:

```
git clone https://github.com/hk-quantum/qday-prize.git
cd qday-prize
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

If you want to run on IBM Quantum hardware, set your IBM Quantum Platform API key and CRN in a `.env` file:

```
API_KEY=<Your API Key>
IBM_CRN=<Your CRN>
```


### How to Run

Specify the bit size to break and the algorithm variant (`compact` or `wide`) as command-line arguments. Examples:

```
python src/main.py 11 compact
```

```
python src/main.py 11 wide
```

The program reads the ECC parameters for the target bit size from `data/curves.json`, constructs the quantum circuit, runs the chosen backend, and writes execution logs. At the end it prints the number of successful private-key recoveries based on 100 measurement shots.

To run on IBM hardware, use the `ibm_main.py` entrypoint with the same arguments (example for 3 bits):

```
python src/ibm_main.py 3 compact
```

```
python src/ibm_main.py 3 wide
```


### Execution Results

From the measurement results of 100 shots, the decrypted private key and the number of successful decryptions are printed in the last line of the log.

The following are example results run on IBM hardware (`ibm_torino`) from [Version 2](https://github.com/hk-quantum/qday-prize/tree/v2.0.0):

|log file|output result|
|---|---|
|`logs/3_ibm.out`|`Success: d=3 count=33`|
|`logs/4_ibm.out`|`Success: d=6 count=20`|


The following results were collected using my simulator `SparseStatevectorSimulator` ([Version 3](https://github.com/hk-quantum/qday-prize/tree/v3.0.0)):

|log file|type|output result|
|---|---|---|
|`logs/3_compsact.out`|`compact`|`Success: d=3 count=75`|
|`logs/3_wide.out`|`wide`|`Success: d=3 count=63`|
|`logs/4_compact.out`|`compact`|`Success: d=6 count=76`|
|`logs/4_wide.out`|`wide`|`Success: d=6 count=81`|
|`logs/6_compact.out`|`compact`|`Success: d=18 count=78`|
|`logs/6_wide.out`|`wide`|`Success: d=18 count=81`|
|`logs/7_compact.out`|`compact`|`Success: d=56 count=86`|
|`logs/7_wide.out`|`wide`|`Success: d=56 count=82`|
|`logs/8_compact.out`|`compact`|`Success: d=103 count=88`|
|`logs/8_wide.out`|`wide`|`Success: d=103 count=89`|
|`logs/9_compact.out`|`compact`|`Success: d=135 count=92`|
|`logs/9_wide.out`|`wide`|`Success: d=135 count=92`|
|`logs/10_compact.out`|`compact`|`Success: d=165 count=91`|
|`logs/11_compact.out`|`compact`|`Success: d=756 count=94`|
|`logs/12_compact.out`|`compact`|`Success: d=1384 count=93`|

For larger bit sizes, only the `compact` variant is shown where the simulator reached its performance limits.
