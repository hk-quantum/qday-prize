# QDays Prize

## Email Address

hk.quantum@icloud.com

## Background

I am a software developer at a Japanese software development company.
To understand how quantum computers work, I have created my own quantum computer simulator and other tools.
For RSA cryptanalysis (factoring a composite number of two primes), I have succeeded in simulating a quantum circuit to break the following bit sizes using my custom quantum circuit simulator.

- Broke 10 bits with Shor's algorithm
- Broke 40 bits with Grover's algorithm

## What I Broke This Time

The maximum bit size for which I was able to actually execute a quantum circuit was 11 bits.

```
--- Bit size 11 ---
Bit size: 11
Prime p: 1051
Curve order (#E): 1093
Subgroup order n: 1093
Cofactor h: 1
Generator point G: (471, 914)
Private key d: 756
Public key Q: (179, 86)
```

## Execution Environment

The quantum circuit was implemented with Qiskit and executed using my high-performance custom quantum computer simulator `SparseStatevectorSimulator` as the backend.

I attempted to execute the circuit on IBM's quantum computer "ibm_torino".  
However, even for a 4-bit ECC curve, the required number of `cz` and `sx` gates exceeded the hardware limit.  
Therefore, the final results are based on simulation.

## Execution Procedure

### Initial Setup

Set up a Python virtual environment and install the necessary libraries.

```
git clone https://github.com/hk-quantum/qday-prize.git
cd qday-prize
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### How to Run

Specify the bit size to be broken as a command-line argument (11 bits in the example).

```
python src/main.py 11
```

The program reads the ECC parameters for the target bit size from `data/curves.json`, builds the quantum circuit, and outputs the execution log of the simulator. Finally, it outputs the number of successful private key decryptions based on 100 shots of measurement.

### Execution Results

From the measurement results of 100 shots, the decrypted private key and the number of successful decryptions are output in the last line of the log.

|log file|output result|
|---|---|
|`logs/4.txt`|`Success: d=6 count=70`|
|`logs/6.txt`|`Success: d=18 count=78`|
|`logs/7.txt`|`Success: d=56 count=95`|
|`logs/8.txt`|`Success: d=103 count=94`|
|`logs/9.txt`|`Success: d=135 count=89`|
|`logs/10.txt`|`Success: d=165 count=90`|
|`logs/11.txt`|`Success: d=756 count=91`|
