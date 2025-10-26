
# Details of the Quantum Circuit

The overall structure of the quantum circuit utilizing Shor's algorithm is shown below.


```
       ┌───┐                         ░ ┌──────┐┌─┐   
   qa: ┤ H ├────────■────────────────░─┤ QFT† ├┤M├───
       ├───┤        │                ░ ├──────┤└╥┘┌─┐
   qb: ┤ H ├────────■────────────────░─┤ QFT† ├─╫─┤M├
       └───┘┌───────┴────────┐┌─┐    ░ └──────┘ ║ └╥┘
   qx: ─────┤0               ├┤M├────░──────────╫──╫─
            │                │└╥┘┌─┐ ░          ║  ║ 
   qy: ─────┤1               ├─╫─┤M├─░──────────╫──╫─
            │  oracle: aG+bQ │ ║ └╥┘ ░          ║  ║ 
anc_0: ─────┤2               ├─╫──╫──░──────────╫──╫─
            │                │ ║  ║  ░          ║  ║ 
anc_1: ─────┤3               ├─╫──╫──░──────────╫──╫─
            └────────────────┘ ║  ║  ░          ║  ║ 
  a: 1/════════════════════════╬══╬═════════════╩══╬═
                               ║  ║             0  ║ 
  b: 1/════════════════════════╬══╬════════════════╩═
                               ║  ║                0 
 xy: 2/════════════════════════╩══╩══════════════════
                               0  1                                  
```

The core of this quantum circuit is the $aG+bQ$ operation. There are two implementations: a "compact" version, which reduces the number of qubits at the cost of deeper circuits, and a "wide" version, which uses more qubits but shallower circuits.

- **Compact version**
       - By uncomputing after each coordinate addition, the number of ancilla qubits is minimized for each ECC addition.
       - Effective for simulators where the number of qubits affects memory usage.  
       ![compcat](aG_add_bQ_compact.png)

- **Wide version**
       - All ECC additions are performed first, then the result is copied and uncomputed at the end.
       - Requires as many qubits as the number of ECC additions, but is effective for real quantum hardware.  
![wide](aG_add_bQ_wide.png)


The set operation $a_iG^{2^i}+b_iQ^{2^i}$ and the ECC addition ($y^2=x^3+ax+b \mod p$) used in the above circuit are described below.

## Quantum Circuit: set $a_iG^{2^i}+b_iQ^{2^i}$

$aG+bQ$ can be decomposed into bitwise additions as $\sum_{i=0}^n (a_iG^{2^i}+b_iQ^{2^i})$.  
For each bit, the addition result is one of four values depending on the combination of $\ket{b_i a_i}$: $\ket{00}O+\ket{01}G^{2^i}+\ket{10}Q^{2^i}+\ket{11}(G^{2^i}+Q^{2^i})$.  
$G^{2^i}, Q^{2^i}, G^{2^i}+Q^{2^i}$ are precomputed classically, and the quantum circuit sets the result to the target register using $b_i, a_i$ as control bits.  
This reduces the number of ECC additions required in the quantum circuit.

```math
\ket{a_i}\ket{b_i}\ket{0} \mapsto \begin{cases} \ket{0}\ket{0}\ket{0}, & \text{if}(a_i=0,b_i=0) \\
\ket{1}\ket{0}\ket{G^{2^i}}, & \text{if}(a_i=1,b_i=0) \\
\ket{0}\ket{1}\ket{Q^{2^i}}, & \text{if}(a_i=0,b_i=1) \\
\ket{1}\ket{1}\ket{G^{2^i}+Q^{2^i}}, & \text{if}(a_i=1,b_i=1) \\
\end{cases}
```

## Quantum Circuit: ECC Addition

ECC addition is a very complex quantum circuit.  
The addition $\ket{xx}+\ket{yy} \mapsto \ket{zz}$ is realized with the following quantum register structure:

```math
\ket{xx}\ket{yy}\ket{zz}\ket{ancilla}
```

Each quantum register is further decomposed as follows:

```math
\begin{align*}
\ket{xx} &= \ket{y_1,x_1} \\
\ket{yy} &= \ket{y_2,x_2} \\
\ket{zz} &= \ket{y_3,x_3} \\
\ket{ancilla} &= \ket{f_1,f_2,f_3,f_4,carry,\lambda,a_1,a_2}
\end{align*}
```

**Details of ancilla:**
- $f_1$: Flag for $(x_1, y_1)=O$
- $f_2$: Flag for $(x_2, y_2)=O$
- $f_3$: Flag for $(x_1, y_1)+(x_2, y_2)=O$
- $f_4$: Flag for $(x_1, y_1)=(x_2, y_2)$ (addition of identical coordinates)
- $carry$: Carry bit used in various calculations
- $\lambda$: Quantum register for $\lambda=\frac{y_2-y_1}{x_2-x_1}$ or $\frac{3x_1^2+a}{2y_1}$ during ECC addition
- $a_1$: Used for $Y$ coordinate calculation in the judgment $(x_1,y_1)+(x_2,y_2)=O$
- $a_2$: Quantum register for modular inverse calculation in $\lambda$ denominator


In Version 3, by efficiently using these registers, the number of quantum gates was greatly reduced.  
The detailed calculation process for ECC addition is as follows:

|Target Register|Control Bits|Quantum Operation|Value After|Remarks|
|---|---|---|---|---|
|$f_1$|$y_1,x_1$|$\ket{0} \mapsto \ket{1}, \text{if}(x_1=0,y_1=0)$|$\begin{cases}\ket{1}, & \text{if}(x_1=0,y_1=0) \\ \ket{0},& \text{otherwise}\end{cases}$|Set flag for $(x_1,y_1)=O$|
|$f_2$|$y_2,x_2$|$\ket{0} \mapsto \ket{1}, \text{if}(x_2=0,y_2=0) $|$\begin{cases}\ket{1}, & \text{if}(x_2=0,y_2=0) \\ \ket{0},& \text{otherwise}\end{cases}$|Set flag for $(x_2,y_2)=O$|
|$x_3$|$f_1, f_2, x_1, x_2$|$\ket{0} \mapsto \begin{cases} \ket{x_2}, & \text{if}(f_1=1) \\ \ket{x_1}, & \text{if}(f_2=1) \end{cases}$|$\begin{cases} \ket{x_2}, & \text{if}(f_1=1) \\ \ket{x_1}, & \text{if}(f_2=1) \\ \ket{0}, & \text{otherwise}\end{cases}$|Set $x_3$ if either $(x_1,y_1)$ or $(x_2,y_2)$ is $O$|
|$y_3$|$f_1, f_2, y_1, y_2$|$\ket{0} \mapsto \begin{cases} \ket{y_2}, & \text{if}(f_1=1) \\ \ket{y_1}, & \text{if}(f_2=1) \end{cases}$|$\begin{cases} \ket{y_2}, & \text{if}(f_1=1) \\ \ket{y_1}, & \text{if}(f_2=1) \\ \ket{0}, & \text{otherwise}\end{cases}$|Set $y_3$ if either $(x_1,y_1)$ or $(x_2,y_2)$ is $O$|
|$x_3$|$f_1,x_2$|$\ket{0} \mapsto \ket{-x_2 \mod p}, \text{if}(f_1 = 0)$|$\begin{cases} \ket{x_2}, & \text{if}(f_1=1) \\ \ket{x_1}, & \text{if}(f_2=1) \\ \ket{-x_2 \mod p}, & \text{otherwise}\end{cases}$|Compute $-x_2$ term in $x_3=\lambda^2-x_1-x_2$|
|$x_2$|$x_1$|$\ket{x_2} \mapsto \ket{x_2-x_1 \mod p}$|$\ket{x_2-x_1 \mod p}$|Compute denominator of $\lambda$|
|$a_1$|$y_1, y_2$|$\ket{0} \mapsto \ket{y_1+y_2-p}$|$\ket{y_1+y_2-p}$|For checking $(x_1,y_1)+(x_2,y_2)=O$|
|$y_2$|$y_1$|$\ket{y_2} \mapsto \ket{y_2-y_1 \mod p}$|$\ket{y_2-y_1 \mod p}$|Compute numerator of $\lambda$|
|$f_3$|$x_2,a_1$|$\ket{0} \mapsto \ket{1}, \text{if}(x_2=0, a_1=0)$|$\begin{cases} \ket{1}, & \text{if}(x_1=x_2, y_1+y_2=p) \\ \ket{0} , & \text{otherwise} \end{cases}$|Set flag for $(x_1,y_1)+(x_2,y_2)=O$|
|$f_4$|$x_2,y_2$|$\ket{0} \mapsto \ket{1}, \text{if}(x_2=0,y_2=0)$|$\begin{cases} \ket{1}, & \text{if}(x_1=x_2,y_1=y_2) \\ \ket{0}, & \text{otherwise} \end{cases}$|Set flag for identical coordinates|
|$x_3$|$f_3, x_1$|$\ket{x_3} \mapsto \ket{x_3+x_1} ,\text{if}(f_3=1)$|$\begin{cases} \ket{x_2}, & \text{if}(f_1=1) \\ \ket{x_1}, & \text{if}(f_2=1) \\ \ket{0}, & \text{if}(f_3=1) \\ \ket{-x_2 \mod p}, & \text{otherwise}\end{cases}$|Set $x_3=0$ if $(x_1,y_1)+(x_2,y_2)=O$|
|$y_2$|$f_4, x_1$|$\ket{0} \mapsto \ket{3 x_1^2+a \mod p} ,\text{if}(f_4=1)$|$\begin{cases} \ket{3 x_1^2+a \mod p} , & \text{if}(f_4=1) \\ \ket{y_2-y_1 \mod p}, & \text{otherwise} \end{cases}$|Numerator of $\lambda$ for identical coordinates|
|$x_2$|$f_4, y_1$|$\ket{0} \mapsto \ket{2 y_1 \mod p}, \text{if}(f_4=1)$|$\begin{cases} \ket{2 y_1 \mod p}, & \text{if}(f_4=1) \\ \ket{x_2-x_1 \mod p}, & \text{otherwise} \end{cases}$|Denominator of $\lambda$ for identical coordinates|
|$a_2$|$x_2$|$\ket{0} \mapsto \ket{x_2^{p-2} \mod p}\ket{\textit{junk}} $|$\begin{cases} \ket{(2y_1)^{-1} \mod p}\ket{\textit{junk}}, & \text{if}(f_4=1) \\ \ket{(x_2-x_1)^{-1} \mod p}\ket{\textit{junk}} , & \text{otherwise} \end{cases}$|Compute modular inverse for denominator of $\lambda$|
|$\lambda$|$y_2, a_2$|$\ket{0} \mapsto \ket{y_2 \cdot a_2 \mod p}$|$ \begin{cases} \ket{(3x_1^2+a) \cdot (2y_1)^{-1} \mod p}, & \text{if}(f_4=1) \\ \ket{(y_2-y_1) \cdot (x_2-x_1)^{-1} \mod p}, & \text{otherwise} \end{cases}$|Compute $\lambda$|
|$x_3$|$f_1, f_2, f_3, x_1$|$\ket{x_3} \mapsto \ket{x_3 + \lambda^2 - x_1 \mod p}, \text{if}(f_1=0,f_2=0,f_3=0) $|$\begin{cases} \ket{x_2}, & \text{if}(f_1=1) \\ \ket{x_1}, & \text{if}(f_2=1) \\ \ket{0}, & \text{if}(f_3=1) \\ \ket{\lambda^2-x_2-x_1 \mod p}, & \text{otherwise}\end{cases}$|Final calculation of $x_3$|
|$x_1$|$x_3$|$\ket{x_1} \mapsto \ket{x_1 - x_3 \mod p}$|$\ket{x_1 - x_3 \mod p}$|Prepare for $y_3$ calculation|
|$y_3$|$f_1,f_2,f_3,x_1,y_1,\lambda$|$\ket{0} \mapsto \ket{\lambda \cdot x_1 - y_1 \mod p}, \text{if}(f_1=0,f_2=0,f_3=0)$|$\begin{cases} \ket{y_2}, & \text{if}(f_1=1) \\ \ket{y_1}, & \text{if}(f_2=1) \\ \ket{0}, & \text{if}(f_3=1) \\ \ket{\lambda (x_1-x_3) - y_1 \mod p}, & \text{otherwise} \end{cases}$|Final calculation of $y_3$|

After obtaining the result, the ECC addition circuit must be uncomputed to return all ancilla to $\ket{0}$.  
In the wide version, uncomputation is performed after all $aG+bQ$ calculations, requiring extra qubits for each addition.  
In the compact version, uncomputation is performed after each addition, allowing efficient reuse of qubits.

Furthermore, in the compact version, the addition is performed as $\ket{xx}\ket{yy} \mapsto \ket{xx+yy}\ket{yy}$, and ancilla bits are returned to $\ket{0}$ as follows:

1. $\ket{xx}\ket{yy}\ket{0}\ket{0}$
2. $\ket{xx}\ket{yy}\ket{0}\ket{xx+yy,\textit{junk}}$: ECC addition
3. $\ket{xx}\ket{yy}\ket{xx+yy}\ket{xx+yy,\textit{junk}}$: Copy result to third register
4. $\ket{xx}\ket{yy}\ket{xx+yy}\ket{0}$: Uncompute ECC addition to return fourth register to $\ket{0}$
5. $\ket{xx+yy}\ket{yy}\ket{xx}\ket{0}$: SWAP $\ket{xx}$ and $\ket{xx+yy}$
6. $\ket{xx+yy}\ket{-yy}\ket{xx}\ket{0}$: Prepare for subtraction by mapping $\ket{yy} \mapsto \ket{-yy}$
7. $\ket{xx+yy}\ket{-yy}\ket{xx}\ket{xx,\textit{junk}}$: ECC addition to compute $\ket{xx}$
8. $\ket{xx+yy}\ket{-yy}\ket{0}\ket{xx,\textit{junk}}$: Return third register to $\ket{0}$
9. $\ket{xx+yy}\ket{-yy}\ket{0}\ket{0}$: Uncompute ECC addition to return fourth register to $\ket{0}$
10. $\ket{xx+yy}\ket{yy}\ket{0}\ket{0}$: Restore $\ket{-yy}$ to original

### Sub-Quantum Circuits


Various sub-quantum circuits are prepared for ECC addition.

- $\ket{x}_n \mapsto \ket{x+1}$
    - Increment circuit for arbitrary bit width. This is the basis for all calculations.
- $\ket{x}_n \mapsto \ket{x+a}$
    - Addition of a constant. Realized by combining increment circuits.
- $\ket{x}_n\ket{y}_m \mapsto \ket{x}\ket{y+x}$
    - Addition between quantum registers. Realized by combining controlled increment circuits.
- $\ket{x}_{n+1} \mapsto \ket{x + a \mod p}$
    - Modular addition of a constant. Implemented so that the carry bit returns to $\ket{0}$.
- $\ket{x}_n\ket{y}_{n+1} \mapsto \ket{x}\ket{y+x \mod p}$
    - Modular addition between quantum registers. Implemented so that the carry bit returns to $\ket{0}$.
- $\ket{x}_{n+1} \mapsto \ket{-x \mod p}$
    - Modular additive inverse of the quantum register itself.
- $\ket{x}_n\ket{0}_n \mapsto \ket{x}\ket{-x \mod p}$
    - Set modular additive inverse to another quantum register.
- $\ket{x}_n\ket{0}_{n+1} \mapsto \ket{x}\ket{x^2 \mod p}$
    - Modular squaring. This is realized by, for all $i,j$, using controlled operations on $\ket{x_i},\ket{x_j}$ and adding the classically precomputed constant $2^{i+j} \mod p$ to the target register. (The same approach is used for subsequent multiplications.)
- $\ket{x}_n\ket{y}_{n+1} \mapsto \ket{x}\ket{y+x^2 \mod p}$
    - Modular squaring and addition at the same time.
- $\ket{x}_n\ket{0}_{n+1} \mapsto \ket{x}\ket{a \cdot x^2 \mod p}$
    - Modular squaring and multiplication by a constant at the same time.
- $\ket{x}_n\ket{y}_{n+1} \mapsto \ket{x}\ket{y+a \cdot x^2 \mod p}$
    - Modular squaring, multiplication by a constant, and addition at the same time.
- $\ket{x}_n\ket{0}_{n+1} \mapsto \ket{x}\ket{a \cdot x \mod p}$
    - Modular multiplication by a constant.
- $\ket{x}_n\ket{y}_{n+1} \mapsto \ket{x}\ket{y+a \cdot x \mod p}$
    - Modular multiplication by a constant and addition at the same time.
- $\ket{x}_n\ket{y}_n\ket{0}_{n+1} \mapsto \ket{x}\ket{y}\ket{x \cdot y \mod n}$
    - Modular multiplication between quantum registers.
- $\ket{x}_n\ket{y}_n\ket{z}_{n+1} \mapsto \ket{x}\ket{y}\ket{z+x \cdot y \mod n}$
    - Modular multiplication and addition at the same time.
- $\ket{x}_n\ket{0}_n\ket{0}_{m \cdot n+1} \mapsto \ket{x}\ket{x^a}\ket{\textit{junk}}$
    - Modular exponentiation by a constant. Realized by combining modular squaring and modular multiplication. Used for modular inverse calculation.

Even for similar calculations, the quantum circuit differs between constant operations and register-to-register operations.
Also, the following two circuits yield the same result when $\ket{y}=\ket{0}$, so only the latter is strictly necessary, but the former can be implemented with fewer quantum gates and is used as appropriate.

- $\ket{x}_n\ket{0}_{n+1} \mapsto \ket{x}\ket{a \cdot x \mod p}$
- $\ket{x}_n\ket{y}_{n+1} \mapsto \ket{x}\ket{y+a \cdot x \mod p}$

## Scale of the Quantum Circuit

The following table shows the number of qubits, quantum gates, and circuit depth for each ECC bit size.  
The number of quantum gates is counted as one for each MCX or SWAP gate, so the actual number after transpilation is much larger, but this is sufficient to represent the scale of the circuit (especially for simulator execution time).

|ecc<br>bits|qbits<br>(compact)|gates<br>(compact)|depth<br>(compact)|qbits<br>(wide)|gates<br>(wide)|depth<br>(wide)|
|--:|--:|--:|--:|--:|--:|--:|
|3|50|5,100|4,563|82|3,388|3,026|
|4|75|19,514|18,720|128|12,992|12,467|
|6|123|165,008|162,005|390|94,214|69,363|
|7|145|479,286|473,012|660|261,272|128,935|
|8|181|1,102,960|1,090,685|979|593,664|251,620|
|9|239|2,510,080|2,497,470|1534|1,338,402|665,831|
|10|245|3,322,938|3,305,080|1725|1,758,806|777,520|
|11|291|6,937,706|6,901,693|2316|3,650,848|1,452,787|
|12|341|9,232,914||3031|4,835,710||

As shown above, the compact version increases the number of qubits gradually, while the wide version increases rapidly.  
However, the wide version allows more parallel execution, so as the ECC bit size increases, the circuit depth increases more slowly compared to the number of gates.

For actual quantum computers, the wide version is more effective, but with the current error rates, it is difficult to execute even 3-bit ECC without errors.
