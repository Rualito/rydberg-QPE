Firstly, consider an Hamiltonian $H^{(1)}$
$$ H^{(1)} = \sum_i H^{1r}_i $$
that flips the metastable $|1\rangle$ state into interacting $|r\rangle$. This can be done through a global field.

When the qubits are on the $|r\rangle$ state, they are allowed to interact via the rydberg interaction $H^{(2)}$, which has the form
$$ H^{(2)} = \sum_{ij} V_{ij} |rr \rangle_{ij} \langle rr| $$

Then, the idea is to flip the control qubit back to $|1\rangle$ and revert the phases of the target qubits. For instance, in a 3 qubit system, we have the Hamiltonians $H_I, H_{II}$ . First, $H^{(1)}$ is applied to enable all-to-all interaction. The all-to-all Hamiltonian is $H_{I}$
$$ H_I=V_{12} |rr\rangle_{12}\langle rr|+V_{13} |rr\rangle_{13}\langle rr|+V_{23} |rr\rangle_{23}\langle rr| $$
which is applied during $\tau_I$, giving the evolution
$$\begin{split} \exp \left(-i\, H_I\,\tau_I \right) = \exp \left\{ -i \,\theta^I_{12} |rr \rangle_{12} \langle rr | \right\} \\ \times \exp \left\{ -i \,\theta^I_{13} |rr \rangle_{13} \langle rr | \right\}\\\times \exp \left\{ -i \,\theta^I_{23} |rr \rangle_{23} \langle rr | \right\} \end{split}  $$
with $\theta_{ij}^I = V_{ij} \tau_I$. Then, if we, for instance, turn off the interaction on qubit 3, the control qubit, by applying $H_3^{1r}$ with a laser field, we get a correlated evolution by the remaining qubits. Then
$$ \exp(-i H_{II} \tau_{II}) = \exp\left\{-i \theta^{II}_{12} |rr\rangle_{12} \langle rr| \right\} $$
Setting $\theta^{II}_{12}=-\theta^{I}_ {12}$, we then invert the phase gained from the interaction between 1 and 2, our target qubits. The final step is to put qubits 1 and 2 back to $|1\rangle$, by applying $H^{1r}_{1}+H^{1r}_{2}$. 


The challenge now became the following: How do we tune $\theta_{ij}^{I, II}$ based on the Rydberg interaction such that we get the desired phases for "Phase In" (more specifically, phases relevant to QFT)
 - Do/Can we use the geometry of the lattice? (yes, but... [[Lattice Geometry Approach]])
 - Can we use different orbitals for the rydberg state? (maybe?)
 -  Multiple rydberg states?
 - Can we orient the orbitals in different directions?
 - Invert the interaction? (this is digital-analog - not great)
 - Tricks from pulse engineering or optics?



$$ \sum_{ \substack{i<j\\i,j\ne c} } $$