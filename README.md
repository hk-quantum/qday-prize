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

The maximum bit size for which I was able to execute a quantum circuit on the simulator was 11 bits.

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

The quantum circuit was implemented with Qiskit and executed using my high-performance custom quantum computer simulator `SparseStatevectorSimulator` as the backend.

On IBM's quantum computer "ibm_torino", I was able to run a 4-bit ECC curve, but due to noise, the results were random and did not achieve the expected accuracy. The execution time for 100 shots was about 35 seconds.

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

Set your IBM Quantum Platform API key and CRN in a `.env` file:

```
API_KEY=<Your API Key>
IBM_CRN=<Your CRN>
```


### How to Run

Specify the bit size to be broken as a command-line argument (11 bits in the example):

```
python src/main.py 11
```

The program reads the ECC parameters for the target bit size from `data/curves.json`, builds the quantum circuit, and outputs the execution log of the simulator. Finally, it outputs the number of successful private key decryptions based on 100 shots of measurement.

To run on an actual IBM quantum computer, use the following command (example for 3 bits):

```
python src/ibm_main.py 3
```


### Execution Results

From the measurement results of 100 shots, the decrypted private key and the number of successful decryptions are output in the last line of the log.

|log file|output result|backend|
|---|---|---|
|`logs/3_ibm.txt`|`Success: d=3 count=33`|ibm_torino|
|`logs/4_ibm.txt`|`Success: d=6 count=20`|ibm_torino|
|`logs/3.txt`|`Success: d=3 count=68`|SparseStatevectorSimulator|
|`logs/4.txt`|`Success: d=6 count=65`|SparseStatevectorSimulator|
|`logs/6.txt`|`Success: d=18 count=83`|SparseStatevectorSimulator|
|`logs/7.txt`|`Success: d=56 count=83`|SparseStatevectorSimulator|
|`logs/8.txt`|`Success: d=103 count=85`|SparseStatevectorSimulator|
|`logs/9.txt`|`Success: d=135 count=84`|SparseStatevectorSimulator|
|`logs/10.txt`|`Success: d=165 count=90`|SparseStatevectorSimulator|
|`logs/11.txt`|`Success: d=756 count=88`|SparseStatevectorSimulator|
