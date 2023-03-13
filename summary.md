# Meeting summaries



## Tuesday, 28th Feb

 - Bernardo and Gabriel shared a few references regarding Quantum Walks and Rydberg Atoms. Slide 2 [here](link)
 - Discussed other problems that might be interesting to tackle for Bernardo's thesis.
    - Grover search with multi-qubit gates. 
    - Optimization algorithms, such as QAOA or VQE, to tackle chemistry problems.
 - Raúl continued his explorations about Quantum Phase Estimation with Rydberg atoms. Slides [here](link)
    - An article was discussed which proposed a logarithmic depth for Quantum Fourier Transform (QFT), using parallel gates. It turns out that the Rydberg platform isn't the most suitable for paralellization, due to unwanted crosstalks. 
    - The QFT circuit can be well approximated using limited size MultiPhase gate (size k < qubits n). Implementing larger MultiPhase gates incurs in a increased error. This error can be minimized with smaller MultiPhase gates (smaller k), at a cost of a small approximation error. A simple noise model for the MultiPhase gate was used to determine the ideal k for a range of n.  
    - The next steps of this work were discussed. More efficient implementations of the MultiPhase gate should be investigated. It also makes sense to run realistic simulations of the process as an argument of its efficiency in current hardware. It could also be interesting to collaborate with experimentalists and showcase the algorithm in real hardware.  



## Tuesday, 7th March

Seminar by Hannes Pichler, from the University of Innsbruck.  
*Title*: From Many-Body Physics to Quantum Optimization with Rydberg Atoms
 - Hannes started by giving an overview of the basic physics of Rydberg atoms.
    - Tweezer arrays are used experimentally to trap and control cold atoms. The loading of cold atoms isn't an ideal process, and a few atoms are lost. Tweezers are used to move the cold atoms into a desired configuration, essentially removing the entropy of the loading process.
    - Single particle Hamiltonians consider the interaction of the valence electron (usually are Alkali atoms) with precisely tuned lasers. There are two parameters, the laser intensity $\Omega$ and the detuning $\Delta$ (this is essentially the frequency difference with the resonant frequency of the Rydberg state).
    - Atoms in the Rydberg state are prone to strong dipole-dipole interactions. These are van der Waals forces, with $V\sim R^{-6}$ dependence. The blockade radius $R_B$ is the point where the laser intensity isn't strong enough to excite an electron to the Rydberg state, given an electron on a Rydberg state within that radius. In essence, it is the point when $V(R_B)=\Omega$.  
    - The Rydberg state is not trapped by a tweezer, but the ground state is. The measurement is done through this phenomena. The atoms that cannot be 'seen' by the laser must have been on the Rydberg state.
 - Spin Liquids, Quantum Kibble-Zurek Mechanism, Dimmer models with Rydberg atoms arrays. .....(request Slides)
 - Quantum Optimization over MIS - Adiabatic Optimization, adiabatic vs variational algorithm. Quasi-Adiabatic algorithm
 - Maximum Independent Set with arbitrary graphs ('upgrade' from unit disk). Crossover gadgets.






