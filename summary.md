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
