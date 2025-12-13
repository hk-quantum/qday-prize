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

The central operation of this circuit is the calculation of the elliptic-curve point combination aG + bQ. Below we explain how this is constructed; as an example we describe the 8-bit case.

If the integers a and b are represented by bit strings

\(a=\ket{a_7a_6\dots a_0},\quad b=\ket{b_7b_6\dots b_0}\),

then

```math
aG+bQ=\ket{a_7}G^{2^7}+\ket{a_6}G^{2^6}+\cdots+\ket{a_0}G+\ket{b_7}Q^{2^7}+\cdots+\ket{b_0}Q
```

That is, we construct the superposition of point additions by adding the contributions of each bit. In the naive decomposition this would require 15 point additions for the 8-bit example.

Because ECC point addition is a relatively heavy quantum subroutine, Version 3 of this work reduced the number of quantum point additions by precomputing small superpositions classically. Concretely, Version 3 computed

```math
aG+bQ=\ket{a_7b_7}f(a_7,b_7)+\ket{a_6b_6}f_2(a_6,b_6)+\cdots+\ket{a_0b_0}f(a_0,b_0)
```

where classically computed maps $f(a_i,b_i) = \ket{a_i}G^{2^i}+\ket{b_i}Q^{2^i}$ (a superposition of up to four coordinates) are loaded into quantum registers. This reduces the number of quantum point additions from 15 to 7.

In Version 4 (this submission) we further reduce the number of quantum point additions by grouping multiple bits and preparing larger precomputed superpositions classically. The decomposition used here is

```math
\begin{aligned}
aG + bQ =\;&
\ket{a_1 a_2 a_3} g_3(a_1,a_2,a_3)
+ \ket{a_4 a_5 a_6} g_3(a_4,a_5,a_6)
+ \ket{a_0 a_7} g_2(a_0,a_7) \\
&\quad
+ \ket{b_1 b_2 b_3} q_3(b_1,b_2,b_3)
+ \ket{b_4 b_5 b_6} q_3(b_4,b_5,b_6)
+ \ket{b_0 b_7} q_2(b_0,b_7)
\end{aligned}
```

Here the classically prepared functions create the following superpositions which are then loaded into quantum registers:

```math
\begin{aligned}
 g_3(a_i,a_j,a_k) &= \ket{a_i}G^{2^i}+\ket{a_j}G^{2^j}+\ket{a_k}G^{2^k} \\
 g_2(a_i,a_j)    &= \ket{a_i}G^{2^i}+\ket{a_j}G^{2^j} \\
 q_3(b_i,b_j,b_k) &= \ket{b_i}Q^{2^i}+\ket{b_j}Q^{2^j}+\ket{b_k}Q^{2^k} \\
 q_2(b_i,b_j)    &= \ket{b_i}Q^{2^i}+\ket{b_j}Q^{2^j}
\end{aligned}
```

Using this approach, the number of quantum point additions is reduced to five.

As an additional micro-optimization, when particular control branches that add complexity (such as doubling the same coordinate or adding the inverse leading to the point at infinity O) are provably unnecessary, we replace the point-addition circuit with a simpler variant. Concretely we use two circuit modes:

```math
\begin{aligned}
add_{normal} &:\; \ket{a_1 a_2 a_3} g_3(a_1,a_2,a_3)
+ \ket{a_4 a_5 a_6} g_3(a_4,a_5,a_6) \\
add_{nodbl} &:\; \ket{a_1a_2\dots a_6}g_6(a_1,a_2,\dots,a_6) + \ket{a_0a_7}g_2(a_0,a_7)
\end{aligned}
```

When creating the superposition for aG, we select between $add_{normal}$ (no doubling or inverse additions occur) and $add_{nodbl}$ (inverse additions can occur but doubling does not). This choice reduces the number of quantum gates; the same logic is applied symmetrically when building bQ.

## Two circuit variants

This implementation provides two variants of the quantum circuit:

- compact: minimizes the number of qubits at the cost of deeper circuits.
  - The circuit uncomputes after each point addition, allowing ancilla reuse.
  - Effective on simulators where memory scales with qubit count.
  - Illustration:
    ![compact](aG_add_bQ_compact.png)

- wide: increases qubit count to reduce circuit depth.
  - All ECC additions are performed (in parallel where possible), then the state is copied and finally uncomputed.
  - Effective for real quantum hardware that benefits from parallelism.
  - Illustration:
    ![wide](aG_add_bQ_wide.png)

With current hardware noise levels neither variant produced the expected measurement outcomes, but we expect that future improvements in quantum hardware will enable correct execution for one of these approaches.

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


In Version 3 and later, by efficiently using the registers described above, we achieved a significant reduction in the number of quantum gates.  
The detailed computation steps are shown below:

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
In the wide variant uncomputation is performed after all $aG+bQ$ calculations, which requires extra qubits proportional to the number of intermediate additions.  In the compact variant each addition is uncomputed immediately, enabling efficient qubit reuse.

The compact variant performs the addition as $\ket{xx}\ket{yy}\ket{0}\ket{0} -> \dots -> \ket{xx+yy}\ket{yy}\ket{0}\ket{0}$ while ensuring ancilla return to $\ket{0}$ by the following sequence:

1. $\ket{xx}\ket{yy}\ket{0}\ket{0}$
2. $\ket{xx}\ket{yy}\ket{0}\ket{xx+yy,\textit{junk}}$ (ECC addition)
3. $\ket{xx}\ket{yy}\ket{xx+yy}\ket{xx+yy,\textit{junk}}$ (copy result to third register)
4. $\ket{xx}\ket{yy}\ket{xx+yy}\ket{0}$ (uncompute ECC addition to return fourth register to $\ket{0}$)
5. $\ket{xx+yy}\ket{yy}\ket{xx}\ket{0}$ (SWAP to move result)
6. $\ket{xx+yy}\ket{-yy}\ket{xx}\ket{0}$ (map $\ket{yy}\mapsto\ket{-yy}$ to prepare subtraction)
7. $\ket{xx+yy}\ket{-yy}\ket{xx}\ket{xx,\textit{junk}}$ (ECC addition to recompute $\ket{xx}$)
8. $\ket{xx+yy}\ket{-yy}\ket{0}\ket{xx,\textit{junk}}$ (clear third register)
9. $\ket{xx+yy}\ket{-yy}\ket{0}\ket{0}$ (uncompute to clear fourth register)
10. $\ket{xx+yy}\ket{yy}\ket{0}\ket{0}$ (restore $\ket{-yy}$ to original)

Note that for additions where doubling or inverse-based addition cannot result in the point at infinity O (\(add_{normal}\), \(add_{nodbl}\)), we omit unnecessary flags and conditional branches to reduce the number of quantum gates.

### Subcircuits

The implementation contains a set of reusable subcircuits used throughout the ECC arithmetic:

- Increment: $\ket{x}_n\mapsto\ket{x+1}$ — arbitrary-width increment.
- Constant addition: $\ket{x}_n\mapsto\ket{x+a}$ — implemented from increment building blocks.
- Register addition: $\ket{x}_n\ket{y}_m\mapsto\ket{x}\ket{y+x}$ — built from controlled increments.
- Modular constant addition: $\ket{x}_{n+1}\mapsto\ket{x+a\bmod p}$ — implemented such that carries are returned to $\ket{0}$.
- Modular register addition: $\ket{x}_n\ket{y}_{n+1}\mapsto\ket{x}\ket{y+x\bmod p}$ — with carry cleared.
- Modular negation: $\ket{x}_{n+1}\mapsto\ket{-x\bmod p}$.
- Set inverse into another register: $\ket{x}_n\ket{0}_n\mapsto\ket{x}\ket{-x\bmod p}$.
- Modular squaring: $\ket{x}_n\ket{0}_{n+1}\mapsto\ket{x}\ket{x^2\bmod p}$ — implemented by controlled adds of classically precomputed constants $2^{i+j}\bmod p$.
- Combined modular squaring and addition/multiplication primitives used to build higher-level routines (multiplication, exponentiation, modular inverse, etc.).

Some circuits are specialized for constant vs register operands; when $\ket{y}=\ket{0}$ the two variants produce the same result, but the constant-operated variant often uses fewer gates and is preferred when applicable.

## Scale of the Quantum Circuit

The table below summarizes the number of qubits, quantum gates and depth for each ECC bit size.  
The number of quantum gates is counted as one for each MCX or SWAP gate; after transpilation the real gate counts on hardware are much larger, but these numbers capture circuit scale for simulator runtime and resource estimation.

|ecc<br>bits|qbits<br>(compact)|gates<br>(compact)|depth<br>(compact)|qbits<br>(wide)|gates<br>(wide)|depth<br>(wide)|
|--:|--:|--:|--:|--:|--:|--:|
|3|47|3,377|3,078|76|1,720|1,542|
|4|71|12,971|12,548|120|6,522|6,283|
|5|93|47,986|46,860|223|24,010|23,412|
|6|117|138,286|136,175|366|62,292|40,528|
|7|138|338,976|335,086|618|164,926|124,419|
|8|163|825,873|817,349|923|387,780|229,635|
|9|230|1,641,514|1,634,177|1462|784,956|468,072|
|10|235|2,287,611|2,276,484|1635|1,115,400|749,091|
|11|280|4,992,898|4,968,256|2206|2,391,440|1,358,478|
|12|329|6,027,756|6,010,187|2899|2,901,368|1,651,805|

As shown, the compact variant increases qubit count slowly while the wide variant grows rapidly.  
However, the wide variant allows more parallel execution, so as the ECC bit size increases, the circuit depth increases more slowly compared to the number of gates.

In practice the wide variant would be favorable for hardware execution, but with the current error rates even a 3-bit ECC instance is difficult to run error-free on available devices.

As a reference, the table below shows the scale of quantum gate counts after transpilation when running the circuits on real hardware (ibm_fez).

|bits|type|qbits|gates|depth|sx gate|cz gate|rz gate|x gate|
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|3|compact|47|1,012,397|545,758|536,619|253,993|220,752|1,021|
|3|wide|76|505,461|270,058|267,887|127,056|109,939|567|
|4|compact|71|4,167,068|2,246,123|2,217,253|1,048,828|897,179|3,794|
|4|wide|120|2,080,814|1,117,623|1,107,430|524,356|447,024|1,990|
|5|compact|93|17,462,887|9,461,983|9,318,314|4,412,259|3,720,215|12,081|

Beyond these entries we could not execute further runs because we exceeded the hardware limits on qubit count or gate count.
