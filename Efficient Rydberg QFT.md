## Recall QFT

With $n=4$ qubits, the QFT circuit looks like this
![[QFT.png]]

The rotations can be grouped in different ways. Either you group them as they are, with each $\widetilde R_i$ acting with the rotation on qubit $i$ and control with the rest: 
![[QFT1-000001.png]]

or move them around, where each qubit $i$ controls the rotations of the qubits with label $< i$
![[QFT2-000001.png]]

There may be an efficient implementation of QFT that scales only with $n$ (number of qubits) by using "Phase Out ($\widetilde S_i$)/In ($\widetilde R_i$)" rydberg gates.

This implementation may exist if we find a way to implement either $\widetilde R_i$ or $\widetilde S_i$ non-sequentially (with all the rotation gates acting at the same time, independently of the number of qubits). 



## $C_kU^m$? 
There is an article discussing the implementation of $C_kU^m$ gates efficiently in Rydberg atoms, might be an interesting approach to start with
[[@liMultipleQubitCkUmLogic2022]]

