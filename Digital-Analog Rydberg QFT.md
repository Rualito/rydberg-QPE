Main reference: [[@martinDigitalanalogQuantumAlgorithm2020]] 

## Digital-Analog QFT

They start by assuming a native all-to-all connectvity 
$$ H_{int} = g\sum_{j<k} Z^{j} Z^{k} $$
For their implementation of the QFT algorithm, they require custom connectivity, so they construct $H_{ZZ}$ by applying $X$ gates on each side. With this scheme they are able to implement
$$ H_{ZZ} = \sum_{k\ne j} g_{jk} Z^{j}Z^{k} $$
The interaction Hamiltonian $H_{int}$ becomes
$$ H_{int} \rightarrow \frac{1}{t_f} \sum_{n<m}^{N} t_{nm} X^n X^m H_{int} X^n X^m  $$
$$ = \frac{1}{t_F} \sum_{j<k}^N \left( \sum_{n<m}^N t_{nm} (-1)^{\delta_{nj}+\delta_{nk}+\delta_{mj}+\delta_{mk}} g \right) Z^i Z^k $$
The term in parenthesis can be compared to $g_{jk}$, and thus we get the expression
$$ g_{jk} = \sum_{n<m}\left[(-1)^{\delta_{nj}+\delta_{nk}+\delta_{mj}+\delta_{mk}}\right] t_{nm} \frac{g}{t_F} $$
This expression can be inverted and thus obtain $t_{nm}$. Therefore a custom coupling Hamiltonian can be obtained from the digital analog protocol. The remainer of the protocol described on the article are just single qubit gates.


## Compatability of Rydberg interactions

The native rydberg interaction is of the form
$$ H_{int} = \sum_{k\ne j} h_{jk} Z^j Z^k$$
where the interaction strength is given by 
$$ h_{jk} \sim \frac{1}{|\vec r_j - \vec r_k|^p} $$with $\vec r_j$ the position of atom $j$ and $p$ depends on the type of interaction, which is tuned by the experiment setup ($p=3$ for dipole-dipole, $p=6$ for Van der Waals). Thus, the interaction Hamiltonian depends on the lattice. The same calculation can be made as above, and we get the following expression
$$ t_F \frac{g_{jk}}{h_{jk}}=\sum_{n<m} M_{jk\,nm}\ t_{nm} $$
$t_{nm}$ can be obtained by vectorizing the indices and inverting $M_{jk\,nm}$. From this expression we can have an ideal of how the analog time scales with the lattice. As the distance between the atoms  $|\vec r_j - \vec r_k|$ increases, the analog time will scale polynomially depending on $p$. The number of atoms $N$ inside a certain radius $R_q$ on a 2D lattice scales with $N\sim R_q^2$, but the analog time scales as $R_q^p$. This means that the time complexity of the QFT algorithm, for a given $N$, is $>O(N^{p/2})$. The analog version already requires $N^2$ steps, so the overall time complexity is $O(N^{2+p/2})$. On a 3D lattice the time complexity would be $O(N^{2+p/3})$.   
