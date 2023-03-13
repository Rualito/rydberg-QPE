
import itertools
import warnings
from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import asdict
from typing import Any, Optional, Union, cast

import matplotlib.pyplot as plt
import numpy as np
import qutip
from numpy.typing import ArrayLike

import pulser.sampler as sampler
from pulser import Pulse, Sequence
from pulser.register import QubitId
from pulser.sampler.samples import _TargetSlot
from pulser.sequence._seq_drawer import draw_sequence
from pulser_simulation.simconfig import SimConfig
from pulser_simulation.simresults import (
    CoherentResults,
    NoisyResults,
    SimulationResults,
)

from pulser.devices._device_datacls import Device  

from pulser_simulation.simulation import Simulation

from generalSequence import GeneralSequence



SUPPORTED_NOISE = {
    "ising": {
        "dephasing",
        "doppler",
        "amplitude",
        "SPAM",
        "depolarizing",
        "eff_noise",
    },
    "XY": {"SPAM"},
    "multi": {
        "dephasing",
        "doppler",
        "amplitude",
        "SPAM",
        "depolarizing",
        "eff_noise",
    }
}
# TODO: How to implement noise in 'multi' mode????

class GeneralSimulation(Simulation):
    # __init__ doesn't look problematic
    # Methods to overload:
    # run()
    # 
    # 
    #  

    def __init__(self, sequence: GeneralSequence, sampling_rate: float = 1, config: Optional[SimConfig] = None, evaluation_times: Union[float, str, ArrayLike] = "Full", with_modulation: bool = False) -> None:
        # super().__init__(sequence, sampling_rate, config, evaluation_times, with_modulation)
        
        if not isinstance(sequence, Sequence):
            raise TypeError(
                "The provided sequence has to be a valid "
                "pulser.Sequence instance."
            )
        if sequence.is_parametrized() or sequence.is_register_mappable():
            raise ValueError(
                "The provided sequence needs to be built to be simulated. Call"
                " `Sequence.build()` with the necessary parameters."
            )
        if not sequence._schedule:
            raise ValueError("The provided sequence has no declared channels.")
        if all(
            sequence._schedule[x][-1].tf == 0
            for x in sequence.declared_channels
        ):
            raise ValueError(
                "No instructions given for the channels in the sequence."
            )
        self._seq = sequence
        self._interaction = 'ising'
        if self._seq._in_xy:
            self._interaction = 'XY'
        elif self._seq._in_multi: 
            self._interaction = 'multi'
        
        self._qdict = self._seq.qubit_info
        self._size = len(self._qdict)
        self._modulated = with_modulation
        if self._modulated and sequence._slm_mask_targets:
            raise NotImplementedError(
                "Simulation of sequences combining an SLM mask and output "
                "modulation is not supported."
            )
        self._tot_duration = self._seq.get_duration(
            include_fall_time=self._modulated
        )
        self.samples_obj = sampler.sample(
            self._seq,
            modulation=self._modulated,
            # The samples are extended by 1 to improve the ODE
            # solver convergence
            extended_duration=self._tot_duration + 1,
        )

        # Type hints for attributes defined outside of __init__
        self.basis_name: str
        self._config: SimConfig
        self.op_matrix: dict[str, qutip.Qobj]
        self.basis: dict[str, qutip.Qobj]
        self.dim: int
        self._eval_times_array: np.ndarray

        if not (0 < sampling_rate <= 1.0):
            raise ValueError(
                "The sampling rate (`sampling_rate` = "
                f"{sampling_rate}) must be greater than 0 and "
                "less than or equal to 1."
            )
        if int(self._tot_duration * sampling_rate) < 4:
            raise ValueError(
                "`sampling_rate` is too small, less than 4 data points."
            )
        self._sampling_rate = sampling_rate
        self._qid_index = {qid: i for i, qid in enumerate(self._qdict)}
        self._collapse_ops: list[qutip.Qobj] = []

        self.sampling_times = self._adapt_to_sampling_rate(
            # Include extra time step for final instruction from samples:
            np.arange(self._tot_duration + 1, dtype=np.double)
            / 1000
        )
        self.evaluation_times = evaluation_times  # type: ignore

        self._bad_atoms: dict[Union[str, int], bool] = {}
        self._doppler_detune: dict[Union[str, int], float] = {}
        # Stores the qutip operators used in building the Hamiltonian
        self.operators: dict[str, defaultdict[str, dict]] = {
            addr: defaultdict(dict) for addr in ["Global", "Local"]
        }
        # Sets the config as well as builds the hamiltonian
        self.set_config(config) if config else self.set_config(SimConfig())
        if hasattr(self._seq, "_measurement"):
            self._meas_basis = self._seq._measurement
        else:
            if self.basis_name in {"digital", "all"}:
                self._meas_basis = "digital"
            else:
                self._meas_basis = self.basis_name
        self.initial_state = "all-ground"

    def set_config(self, cfg: SimConfig) -> None:
        """Sets current config to cfg and updates simulation parameters.

        Args:
            cfg: New configuration.
        """
        if not isinstance(cfg, SimConfig):
            raise ValueError(f"Object {cfg} is not a valid `SimConfig`.")
        not_supported = set(cfg.noise) - SUPPORTED_NOISE[self._interaction]
        if not_supported:
            raise NotImplementedError(
                f"Interaction mode '{self._interaction}' does not support "
                f"simulation of noise types: {', '.join(not_supported)}."
            )
        prev_config = self.config if hasattr(self, "_config") else SimConfig()
        self._config = cfg
        if not ("SPAM" in self.config.noise and self.config.eta > 0):
            self._bad_atoms = {qid: False for qid in self._qid_index}
        if "doppler" not in self.config.noise:
            self._doppler_detune = {qid: 0.0 for qid in self._qid_index}
        # Noise, samples and Hamiltonian update routine
        self._construct_hamiltonian()

        kraus_ops = []
        if "dephasing" in self.config.noise:
            if self.basis_name == "digital" or self.basis_name == "all":
                # Go back to previous config
                self.set_config(prev_config)
                raise NotImplementedError(
                    "Cannot include dephasing noise in digital- or all-basis."
                )
            # Probability of phase (Z) flip:
            # First order in prob
            prob = self.config.dephasing_prob / 2
            n = self._size
            if prob > 0.1 and n > 1:
                warnings.warn(
                    "The dephasing model is a first-order approximation in the"
                    f" dephasing probability. p = {2*prob} is too large for "
                    "realistic results.",
                    stacklevel=2,
                )
            k = np.sqrt(prob * (1 - prob) ** (n - 1))

            self._collapse_ops += [
                np.sqrt((1 - prob) ** n)
                * qutip.tensor([self.op_matrix["I"] for _ in range(n)])
            ]
            kraus_ops.append(k * qutip.sigmaz())

        if "depolarizing" in self.config.noise:
            if self.basis_name == "digital" or self.basis_name == "all":
                # Go back to previous config
                self.set_config(prev_config)
                raise NotImplementedError(
                    "Cannot include depolarizing "
                    + "noise in digital- or all-basis."
                )
            # Probability of error occurrence

            prob = self.config.depolarizing_prob / 4
            n = self._size
            if prob > 0.1 and n > 1:
                warnings.warn(
                    "The depolarizing model is a first-order approximation"
                    f" in the depolarizing probability. p = {4*prob}"
                    " is too large for realistic results.",
                    stacklevel=2,
                )

            k = np.sqrt((prob) * (1 - 3 * prob) ** (n - 1))
            self._collapse_ops += [
                np.sqrt((1 - 3 * prob) ** n)
                * qutip.tensor([self.op_matrix["I"] for _ in range(n)])
            ]
            kraus_ops.append(k * qutip.sigmax())
            kraus_ops.append(k * qutip.sigmay())
            kraus_ops.append(k * qutip.sigmaz())

        if "eff_noise" in self.config.noise:
            if self.basis_name == "digital" or self.basis_name == "all":
                # Go back to previous config
                self.set_config(prev_config)
                raise NotImplementedError(
                    "Cannot include general "
                    + "noise in digital- or all-basis."
                )
            # Probability distribution of error occurences
            n = self._size
            m = len(self.config.eff_noise_opers)
            if n > 1:
                for i in range(1, m):
                    prob_i = self.config.eff_noise_probs[i]
                    if prob_i > 0.1:
                        warnings.warn(
                            "The effective noise model is a first-order"
                            " approximation in the noise probability."
                            f"p={prob_i} is large for realistic results.",
                            stacklevel=2,
                        )
                        break
            # Deriving Kraus operators
            prob_id = self.config.eff_noise_probs[0]
            self._collapse_ops += [
                np.sqrt(prob_id**n)
                * qutip.tensor([self.op_matrix["I"] for _ in range(n)])
            ]
            for i in range(1, m):
                k = np.sqrt(
                    self.config.eff_noise_probs[i] * prob_id ** (n - 1)
                )
                k_op = k * self.config.eff_noise_opers[i]
                kraus_ops.append(k_op)

        # Building collapse operators
        for operator in kraus_ops:
            self._collapse_ops += [
                self.build_operator([(operator, [qid])])
                for qid in self._qid_index
            ]

    def add_config(self, config: SimConfig) -> None:
        """Updates the current configuration with parameters of another one.

        Mostly useful when dealing with multiple noise types in different
        configurations and wanting to merge these configurations together.
        Adds simulation parameters to noises that weren't available in the
        former SimConfig. Noises specified in both SimConfigs will keep
        former noise parameters.

        Args:
            config: SimConfig to retrieve parameters from.
        """
        if not isinstance(config, SimConfig):
            raise ValueError(f"Object {config} is not a valid `SimConfig`")

        not_supported = set(config.noise) - SUPPORTED_NOISE[self._interaction]
        if not_supported:
            raise NotImplementedError(
                f"Interaction mode '{self._interaction}' does not support "
                f"simulation of noise types: {', '.join(not_supported)}."
            )

        old_noise_set = set(self.config.noise)
        new_noise_set = old_noise_set.union(config.noise)
        diff_noise_set = new_noise_set - old_noise_set
        # Create temporary param_dict to add noise parameters:
        param_dict: dict[str, Any] = asdict(self._config)
        # remove redundant `spam_dict`:
        del param_dict["spam_dict"]
        # `doppler_sigma` will be recalculated from temperature if needed:
        del param_dict["doppler_sigma"]
        # Begin populating with added noise parameters:
        param_dict["noise"] = tuple(new_noise_set)
        if "SPAM" in diff_noise_set:
            param_dict["eta"] = config.eta
            param_dict["epsilon"] = config.epsilon
            param_dict["epsilon_prime"] = config.epsilon_prime
        if "doppler" in diff_noise_set:
            param_dict["temperature"] = config.temperature
        if "amplitude" in diff_noise_set:
            param_dict["laser_waist"] = config.laser_waist
        if "dephasing" in diff_noise_set:
            param_dict["dephasing_prob"] = config.dephasing_prob
        if "depolarizing" in diff_noise_set:
            param_dict["depolarizing_prob"] = config.depolarizing_prob
        if "eff_noise" in diff_noise_set:
            param_dict["eff_noise_opers"] = config.eff_noise_opers
            param_dict["eff_noise_probs"] = config.eff_noise_probs
        param_dict["temperature"] *= 1.0e6
        # update runs:
        param_dict["runs"] = config.runs
        param_dict["samples_per_run"] = config.samples_per_run

        # set config with the new parameters:
        self.set_config(SimConfig(**param_dict))


    def _build_basis_and_op_matrices(self) -> None:
        self.basis_ids = {}
        if self._interaction == "multi":
            # Digital, ground-rydberg and XY
            self.basis_name = 'multi'
            self.dim = 2 + len(self._seq.device.rydberg_states)
            basis = ['g', 'h'] + [*self._seq.device.state_labels]
            # digital basis + ground-rydbergs + ry-ry 
            projectors = [('g', 'g'), ('h','g')] + [ ('g',ry) for ry in self._seq.device.state_labels ] + [(r1,r2) for r1 in self._seq.device.state_labels for r2 in self._seq.device.state_labels]
            self.basis_ids = { str(('g', 'h')): ['sigma_gg', 'sigma_gh']}
            for proj in projectors:
                s1, s2 = proj
                if s1.__contains__(s2):
                    continue
                self.basis_ids[str( proj )] = ['sigma_'+str(proj), 'sigma_'+str((s2, s2))]
        elif self._interaction == "XY":
            self.basis_name = "XY"
            self.dim = 2
            basis = ["u", "d"]
            projectors = ["uu", "du", "ud", "dd"]
        else:
            used_bases = self.samples_obj.used_bases()
            if "digital" not in used_bases:
                self.basis_name = "ground-rydberg"
                self.dim = 2
                basis = ["r", "g"]
                projectors = ["gr", "rr", "gg"]
            elif "ground-rydberg" not in used_bases:
                self.basis_name = "digital"
                self.dim = 2
                basis = ["g", "h"]
                projectors = ["hg", "hh", "gg"]
            else:
                self.basis_name = "all"  # All three states
                self.dim = 3
                basis = ["r", "g", "h"]
                projectors = ["gr", "hg", "rr", "gg", "hh"]
        
        self.basis = {b: qutip.basis(self.dim, i) for i, b in enumerate(basis)}
        self.op_matrix = {"I": qutip.qeye(self.dim)}

        for proj in projectors:
            self.op_matrix["sigma_" + str(proj)] = (
                self.basis[proj[0]] * self.basis[proj[1]].dag()
            )

    def _construct_hamiltonian(self, update: bool = True) -> None:
        """Constructs the hamiltonian from the Sequence.

        Also builds qutip.Qobjs related to the Sequence if not built already,
        and refreshes potential noise parameters by drawing new at random.

        Args:
            update: Whether to update the noise parameters.
        """
        if update:
            self._update_noise()
        self._extract_samples()
        if not hasattr(self, "basis_name"):
            self._build_basis_and_op_matrices()

        def make_vdw_term(q1: QubitId, q2: QubitId, ry1:str, ry2:str) -> qutip.Qobj:
            # TODO: vdw terms between rr and RR? Device needs two interaction terms
            """Construct the Van der Waals interaction Term.

            For each pair of qubits, calculate the distance between them,
            then assign the local operator "sigma_rr" at each pair.
            The units are given so that the coefficient includes a
            1/hbar factor.
            """
            dist = np.linalg.norm(self._qdict[q1] - self._qdict[q2])
            op = self.build_operator([("sigma_"+str( (ry1,ry2)), [q1, q2]), ("sigma_"+str((ry2,ry1)), [q2, q1])])
            if (ry1, ry2) in self._seq.device.c6_dict:
                U = 0.5 * self._seq.device.c6_dict[(ry1, ry2)]/ dist**6
                return U * op
            
            return 0*op

        def make_xy_term(q1: QubitId, q2: QubitId, ry1:str, ry2:str) -> qutip.Qobj:
            # TODO: AsymSequence has to define a magnetic field
            """Construct the XY interaction Term.

            For each pair of qubits, calculate the distance between them,
            then assign the local operator "sigma_du * sigma_ud" at each pair.
            The units are given so that the coefficient
            includes a 1/hbar factor.
            """
            dist = np.linalg.norm(self._qdict[q1] - self._qdict[q2])
            coords_dim = len(self._qdict[q1])
            mag_norm = np.linalg.norm(self._seq.magnetic_field[:coords_dim])
            if mag_norm < 1e-8:
                cosine = 0.0
            else:
                cosine = np.dot(
                    (self._qdict[q1] - self._qdict[q2]),
                    self._seq.magnetic_field[:coords_dim],
                ) / (dist * mag_norm)
            U = (
                0.5
                * cast(float, self._seq._device.interaction_coeff_xy)
                * (1 - 3 * cosine**2)
                / dist**3
            )
            if self._interaction == 'multi':
                if (ry1, ry2) in self._seq.device.c3_dict:
                    U = (
                        0.5
                        * cast(float, self._seq.device.c3_dict[(ry1, ry2)])
                        * (1 - 3 * cosine**2)
                        / dist**3
                    )
                    return U * self.build_operator([("sigma_"+str( (ry1,ry2)), [q1, q2]), ("sigma_"+str((ry2,ry1)), [q2, q1])])
            
            if self._interaction == 'XY':
                return U * self.build_operator(
                    [("sigma_du", [q1]), ("sigma_ud", [q2])]
                )
            return 0

        def make_interaction_term(masked: bool = False) -> qutip.Qobj:
            if masked:
                # Calculate the total number of good, unmasked qubits
                effective_size = self._size - sum(self._bad_atoms.values())
                for q in self._seq._slm_mask_targets:
                    if not self._bad_atoms[q]:
                        effective_size -= 1
                if effective_size < 2:
                    return 0 * self.build_operator([("I", "global")])

            # make interaction term
            dipole_interaction = cast(qutip.Qobj, 0)
            for q1, q2 in itertools.combinations(self._qdict.keys(), r=2):
                if (
                    self._bad_atoms[q1]
                    or self._bad_atoms[q2]
                    or (
                        masked
                        and (
                            q1 in self._seq._slm_mask_targets
                            or q2 in self._seq._slm_mask_targets
                        )
                    )
                ):
                    continue

                if self._interaction == "XY":
                    dipole_interaction += make_xy_term(q1, q2)
                elif self._interaction == 'ising':
                    dipole_interaction += make_vdw_term(q1, q2)
                elif self._interaction == 'multi':

                    # Go through all possible state combinations
                    for ry1 in self._seq.device.state_labels:
                        for ry2 in self._seq.device.state_labels:
                            dist = np.linalg.norm(self._qdict[q1] - self._qdict[q2])

                            # if radius less than LeRoy radius, then assume XY interaction
                            if dist < self._seq.device.rvdw_dict[(ry1, ry2)]: 
                                dipole_interaction += make_xy_term(q1, q2, ry1, ry2)
                            else:
                                # otherwise assume van der Waals
                                dipole_interaction += make_vdw_term(q1, q2, ry1, ry2)

                
            return dipole_interaction

        def build_coeffs_ops(basis: str, addr: str) -> list[list]:
            """Build coefficients and operators for the hamiltonian QobjEvo."""
            samples = self.samples[addr][basis]
            operators = self.operators[addr][basis]
            # Choose operator names according to addressing:
            # TODO: change op_ids according to basis. All the possible basis come from device/channels. 
            if basis == "ground-rydberg":
                op_ids = ["sigma_gr", "sigma_rr"]
            elif basis == "digital":
                op_ids = ["sigma_hg", "sigma_gg"]
            elif basis == "XY":
                op_ids = ["sigma_du", "sigma_dd"]
            elif basis == 'ground-RYDBERG':
                op_ids = ["sigma_gR", "sigma_RR"]
            elif self._interaction == 'multi':
                op_ids = self.basis_ids[basis]
            # if self._interaction == 'asym':
            #     op_ids += "sigma_rR", "sigma_Rr"
            # elif basis == 'asym':
            #     op_ids = []

            terms = []
            if addr == "Global":
                coeffs = [
                    0.5 * samples["amp"] * np.exp(-1j * samples["phase"]),
                    -0.5 * samples["det"],
                ]
                for op_id, coeff in zip(op_ids, coeffs):
                    if np.any(coeff != 0):
                        # Build once global operators as they are needed
                        if op_id not in operators:
                            operators[op_id] = self.build_operator(
                                [(op_id, "global")]
                            )
                        terms.append(
                            [
                                operators[op_id],
                                self._adapt_to_sampling_rate(coeff),
                            ]
                        )
            elif addr == "Local":
                for q_id, samples_q in samples.items():
                    if q_id not in operators:
                        operators[q_id] = {}
                    coeffs = [
                        0.5
                        * samples_q["amp"]
                        * np.exp(-1j * samples_q["phase"]),
                        -0.5 * samples_q["det"],
                    ]
                    for coeff, op_id in zip(coeffs, op_ids):
                        if np.any(coeff != 0):
                            if op_id not in operators[q_id]:
                                operators[q_id][op_id] = self.build_operator(
                                    [(op_id, [q_id])]
                                )
                            terms.append(
                                [
                                    operators[q_id][op_id],
                                    self._adapt_to_sampling_rate(coeff),
                                ]
                            )
            self.operators[addr][basis] = operators
            return terms

        qobj_list = []
        # Time independent term:
        effective_size = self._size - sum(self._bad_atoms.values())
        if self.basis_name != "digital" and effective_size > 1:
            # Build time-dependent or time-independent interaction term based
            # on whether an SLM mask was defined or not
            if self._seq._slm_mask_time:
                # Build an array of binary coefficients for the interaction
                # term of unmasked qubits
                coeff = np.ones(self._tot_duration)
                coeff[0 : self._seq._slm_mask_time[1]] = 0
                # Build the interaction term for unmasked qubits
                qobj_list = [
                    [
                        make_interaction_term(),
                        self._adapt_to_sampling_rate(coeff),
                    ]
                ]
                # Build the interaction term for masked qubits
                qobj_list += [
                    [
                        make_interaction_term(masked=True),
                        self._adapt_to_sampling_rate(
                            np.logical_not(coeff).astype(int)
                        ),
                    ]
                ]
            else:
                qobj_list = [make_interaction_term()]
        # print(self.samples)
        # Time dependent terms:
        for addr in self.samples:
            for basis in self.samples[addr]:
                if self.samples[addr][basis]:
                    qobj_list += cast(list, build_coeffs_ops(basis, addr))

        if not qobj_list:  # If qobj_list ends up empty
            qobj_list = [0 * self.build_operator([("I", "global")])]

        ham = qutip.QobjEvo(qobj_list, tlist=self.sampling_times)
        ham = ham + ham.dag()
        ham.compress()
        self._hamiltonian = ham


