# Details of the Quantum Circuit


## Details of the Oracle

In Version 2, the number of quantum gates required to realize $aG+bQ$ was significantly reduced. The differences between Version 1 and Version 2 are shown below.

In Version 1, the classically computed $G^1, G^2, G^4$ and $Q^1, Q^2, Q^4$ were sequentially added using the corresponding bits of the quantum registers $\ket{a_2 a_1 a_0}$ and $\ket{b_2 b_1 b_0}$ as control bits.

```
 qa_0: ────■─────────────────────────────────────────────────
           │                                                 
 qa_1: ────┼────────■────────────────────────────────────────
           │        │                                        
 qa_2: ────┼────────┼────────■───────────────────────────────
           │        │        │                               
 qb_0: ────┼────────┼────────┼────────■──────────────────────
           │        │        │        │                     
 qb_1: ────┼────────┼────────┼────────┼────────■─────────────
           │        │        │        │        │            
 qb_2: ────┼────────┼────────┼────────┼────────┼────────■────
       ┌───┴───┐┌───┴───┐┌───┴───┐┌───┴───┐┌───┴───┐┌───┴───┐
   qx: ┤0      ├┤0      ├┤0      ├┤0      ├┤0      ├┤0      ├
       │       ││       ││       ││       ││       ││       │
   qy: ┤1      ├┤1      ├┤1      ├┤1      ├┤1      ├┤1      ├
       │       ││       ││       ││       ││       ││       │
anc_0: ┤2 +G^1 ├┤2 +G^2 ├┤2 +G^4 ├┤2 +Q^1 ├┤2 +Q^2 ├┤2 +Q^4 ├
       │       ││       ││       ││       ││       ││       │
anc_1: ┤3      ├┤3      ├┤3      ├┤3      ├┤3      ├┤3      ├
       │       ││       ││       ││       ││       ││       │
anc_2: ┤4      ├┤4      ├┤4      ├┤4      ├┤4      ├┤4      ├
       └───────┘└───────┘└───────┘└───────┘└───────┘└───────┘
```

On the other hand, in Version 2, $a_i G^{2^i} + b_i Q^{2^i}$ is classically computed for each bit and embedded into the quantum register.

```math
\ket{a_i}\ket{b_i}=\ket{1}\ket{0}G^{2^i}+\ket{0}\ket{1}Q^{2^i}+\ket{1}\ket{1}(G^{2^i}+Q^{2^i})
```

Since there are only three possible combinations of bits, embedding into the quantum register can be done simply.

By using this embedded quantum register for addition, the number of ECC additions, which are complex quantum circuits, was greatly reduced.  
(For example, in the 3-bit case, the number of additions was reduced from 6 to 2)

```
 qa_0: ───────■──────────────────────────────────────────────────────────────────────────────────────────
              │                                                                                          
 qa_1: ───────┼─────────────────■────────────────────────────────────────────────────────────────────────
              │                 │                                                                        
 qa_2: ───────┼─────────────────┼───────────────────■────────────────────────────────────────────────────
              │                 │                   │                                                    
 qb_0: ───────■─────────────────┼───────────────────┼────────────────────────────────────────────────────
              │                 │                   │                                                    
 qb_1: ───────┼─────────────────■───────────────────┼────────────────────────────────────────────────────
              │                 │                   │                                                    
 qb_2: ───────┼─────────────────┼───────────────────■────────────────────────────────────────────────────
       ┌──────┴───────┐         │                   │          ┌───────────────────┐┌───────────────────┐
   qx: ┤0             ├─────────┼───────────────────┼──────────┤0                  ├┤0                  ├
       │  a_0*G+b_0*Q │         │                   │          │                   ││                   │
   qy: ┤1             ├─────────┼───────────────────┼──────────┤1                  ├┤1                  ├
       └──────────────┘         │                   │          │                   ││                   │
anc_0: ─────────────────────────┼───────────────────┼──────────┤2 +a_1*G^2+b_1*Q^2 ├┤2 +a_2*G^4+b_2*Q^4 ├
                                │                   │          │                   ││                   │
anc_1: ─────────────────────────┼───────────────────┼──────────┤3                  ├┤3                  ├
                                │                   │          │                   ││                   │
anc_2: ─────────────────────────┼───────────────────┼──────────┤4                  ├┤4                  ├
                       ┌────────┴─────────┐         │          └─────────┬─────────┘└─────────┬─────────┘
  qx2: ────────────────┤0                 ├─────────┼────────────────────■────────────────────┼──────────
                       │  a_1*G^2+b_1*Q^2 │         │                    │                    │          
  qy2: ────────────────┤1                 ├─────────┼────────────────────■────────────────────┼──────────
                       └──────────────────┘┌────────┴─────────┐                               │          
  qx4: ────────────────────────────────────┤0                 ├───────────────────────────────■──────────
                                           │  a_2*G^4+b_2*Q^4 │                               │          
  qy4: ────────────────────────────────────┤1                 ├───────────────────────────────■──────────
                                           └──────────────────┘                                          
```

This method not only reduces the number of quantum gates, but also allows the initial embedding into the quantum register to be performed in parallel, thus reducing the circuit depth.

## Addition of Coordinates

The addition of complex coordinates is realized by the following quantum circuit.

```
  qx1: ─────────■───────────X───────────────────────■────────────────────
                │           │                       │                    
  qy1: ─────────■───────────┼──X────────────────────■────────────────────
                │           │  │                    │                    
  qx2: ─────────■───────────┼──┼────────────────────■────────────────────
                │           │  │ ┌──────┐           │            ┌──────┐
  qy2: ─────────■───────────┼──┼─┤ -qy2 ├───────────■────────────┤ -qy2 ├
       ┌────────┴─────────┐ │  │ └──────┘┌──────────┴───────────┐└──────┘
   ox: ┤0                 ├─X──┼─────────┤0                     ├────────
       │                  │    │         │                      │        
   oy: ┤1                 ├────X─────────┤1                     ├────────
       │                  │              │                      │        
anc_0: ┤2 (x1,y1)+(x2,y2) ├──────────────┤2 inv (x1,y1)+(x2,y2) ├────────
       │                  │              │                      │        
anc_1: ┤3                 ├──────────────┤3                     ├────────
       │                  │              │                      │        
anc_2: ┤4                 ├──────────────┤4                     ├────────
       └──────────────────┘              └──────────────────────┘        
```

- The result of adding the inputs $(qx1, qy1)$ and $(qx2, qy2)$ is stored in $(ox, oy)$
- $(qx1, qy1)$ and $(ox, oy)$ are swapped
- For subtraction, $qy2$ is inverted as $-qy2 = p - qy2$
- By running the inverse circuit of $(qx1, qy1)+(qx2, qy2)$ in this state, all of $ox, oy, anc$ are returned to $\ket{0}$

The following ancilla bits are prepared for this addition circuit:

- Four 1-bit flags indicating the following states:
  - $(qx1, qy1)=O$
  - $(qx2, qy2)=O$
  - $(qx1, qy1)=-(qx2, qy2)$
  - $(qx, qy)=(qx2, qy2)$
- $dx:=(qx1-qx2) \mod{p}$
- $dy:=(qy1-qy2) \mod{p}$
- $dx^{-1}:=dx^{p-2} \mod{p}$
  - Additional quantum bits are prepared to hold the powers $dx^2,dx^4,dx^8...$ up to the most significant bit of $p-2$
  - Temporary quantum registers are also prepared for calculation
- $\lambda:=dx^{-1} \cdot dy \mod{p}$
- $\lambda \cdot ox$ (for the calculation of $oy$)
- 1-bit carry bit (used for addr mod in various calculations)
