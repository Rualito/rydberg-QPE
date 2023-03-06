Consider that we have $n=8$ qubits. The QFT circuit, using MultiPhase, is the following:
![[QFT-Efficient-exact-000001.png]]
After applying swaps, this is the resulting state
$$ | \psi_x \rangle = \frac{1}{\sqrt {2^n}}\ \left(\,|0\rangle + e^{2\pi i\ (0.x_0)}|1\rangle\,\right) \cdots \left(\,|0\rangle + e^{2\pi i\ (0.x_{n-1}...x_{0})}|1\rangle\,\right) $$
which can be rewritten as
$$|\psi_{x} \rangle = | \mu_{0.x_0} \rangle | \mu_{0.x_1x_0} \rangle \cdots | \mu_{x_{n-1}...x_0} \rangle$$
where
$$ | \mu_\theta \rangle = \frac{1}{\sqrt 2} \left( |0\rangle + e^{2 \pi i\ \theta} |1\rangle \right) $$
The first stage of the circuit above (Hadamard + MultiPhase) computes the state $|\mu_{0.x_1...x_8} \rangle_1$. Although, notice how it requires that one control qubit to have an effect on 7 target qubits. This is likely to not be realistic. Thus, consider the following circuit, corresponding to an approximation of the QFT.
![[QFT-Efficient-approx-000001.png]]
Now we use MultiPhase gates with at most 3 target qubits, with an approximation cost.

For instance, now the first stage computes $|\mu_{0.x_1...x_4}\rangle$. Thus we have the following error
$$ ||\ |\mu_{0.x_1...x_8}\rangle -|\mu_{0.x_1...x_8}\rangle\  || $$

