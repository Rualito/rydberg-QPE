from sys import version_info

from pulser.devices._device_datacls import Device, BaseDevice  
from dataclasses import dataclass, field, fields
from pulser.json.utils import obj_to_dict

if version_info[:2] >= (3, 8):  # pragma: no cover
    from typing import Literal, get_args
else:  # pragma: no cover
    try:
        from typing_extensions import Literal, get_args  # type: ignore
    except ImportError:
        raise ImportError(
            "Using pulser with Python version 3.7 requires the"
            " `typing_extensions` module. Install it by running"
            " `pip install typing-extensions`."
        )

from generalChannel import GeneralChannel

import arc
DIMENSIONS = Literal[2, 3]

from typing import Any, cast
from collections import Counter

import numpy as np

ratio_ARC_TO_PASQAL = 2*np.pi * 10**3 

from pulser.devices._devices import Chadoq2


def getC3fromstates(atom1, state1_s, state1_f, atom2, state2_s, state2_f):
    '''
    Estimate the short range coeficient C_3.
    atom1,2: atoms that are interacting
    state1,2_s: states of corresponding atoms where the interaction is starting
    state1,2_f: states of corresponding atoms where the interaction ends
    q: specifies transition that the driving field couples to, +1, 0 or -1 corresponding to driving
    sigma_+, pi, sigma_- transitions respectively.   
    Returns C_3 in units of MHz (mu m)^3 (adjusted to PasQal units)
    '''
    # coupling od 59 D_{3/2} m_j = 3/2 -> 51 P_{1/2} m_j = 1/2

    dpDME = atom1.getDipoleMatrixElement(*state1_s, *state1_f, state1_f[3]-state1_s[3])
    # coupling od 59 D_{3/2} m_j = 3/2 -> 57 F_{5/2} m_j = 5/2
    dfDME = atom2.getDipoleMatrixElement(*state2_s, *state2_f, state2_f[3]-state2_s[3])
    c3  = 1/(4.0*np.pi*arc.epsilon_0)*dpDME*dfDME*arc.C_e**2*\
                    (arc.physical_constants["Bohr radius"][0])**2
    return (abs(c3)/arc.C_h*1.e9) * ratio_ARC_TO_PASQAL
    # print("C_3 = %.3f GHz (mu m)^3 " % (abs(c3)/arc.C_h*1.e9  ))


def getC6fromstates(atom1, state1, atom2, state2, theta=np.pi/2):
    '''
    Estimate the long range coeficient C_6.
    atom1,2: atoms that are interacting
    state1,2: states of corresponding atoms that are interacting 
    Returns C_6 in units of MHz (mu m)^6 (adjusted to PasQal units)
    '''
    # parameters: atom, n1, l1, j1 , n2, l2 , j2 , mj1, mj2 
    calc = arc.PairStateInteractions(atom1, *state1[:3], *state2[:3], state1[3], state2[3], atom2=atom2) 
    # parameters: theta, phi, deltan, deltaE 
    c6 = calc.getC6perturbatively(theta, 0, 5, 25.e9) 
    # print("C_6 = %.0f GHz (mu m)^6" % c6)
    return c6 * ratio_ARC_TO_PASQAL

def getRvdw(atom1, state1, atom2, state2):
    '''
    Compute the LeRoy radius for the particular atom-state combination.
    '''
    calc = arc.PairStateInteractions(atom1, *state1[:3], *state2[:3], state1[3], state2[3], atom2=atom2) 
    return calc.getLeRoyRadius()



@dataclass(frozen=True, repr=False)  # type: ignore[misc]
class GeneralDevice(BaseDevice):
    # Define C6 for arbitrary atoms
    # Compatible with the general framework
    name: str
    dimensions: DIMENSIONS
    min_atom_distance: float = Chadoq2.min_atom_distance
    max_atom_num: int = Chadoq2.max_atom_num
    max_radial_distance: int | None = Chadoq2.max_radial_distance
    interaction_coeff_xy: float | None = None
    supports_slm_mask: bool = False
    alkali_atom: arc.AlkaliAtom = arc.Rubidium()
    rydberg_states: tuple[tuple[int]] = ( (60,0,0.5,0.5), (60,0,0.5,-0.5)) 
    state_labels: tuple[str]= ('r1', 'r2')
    rydberg_level: int = 0
    max_layout_filling: float = 0.5
    reusable_channels: bool = field(default=False, init=False)
    channel_ids: tuple[str, ...] | None = None
    channel_objects: tuple[GeneralChannel, ...] = field(default_factory=tuple)
    c3_dict: dict = field(default_factory=dict)
    c6_dict: dict = field(default_factory=dict)
    rvdw_dict: dict = field(default_factory=dict)
    

    def __post_init__(self) -> None:
        def type_check(
            param: str, type_: type, value_override: Any = None
        ) -> None:
            value = (
                getattr(self, param)
                if value_override is None
                else value_override
            )
            if not isinstance(value, type_):
                raise TypeError(
                    f"{param} must be of type '{type_.__name__}', "
                    f"not '{type(value).__name__}'."
                )

        type_check("name", str)
        if self.dimensions not in get_args(DIMENSIONS):
            raise ValueError(
                f"'dimensions' must be one of {get_args(DIMENSIONS)}, "
                f"not {self.dimensions}."
            )

        self._validate_rydberg_states(self.rydberg_states, self.state_labels)
        

        for param in (
            "min_atom_distance",
            "max_atom_num",
            "max_radial_distance",
        ):
            value = getattr(self, param)
            if param in self._optional_parameters:
                prelude = "When defined, "
                is_none = value is None
            elif value is None:
                raise TypeError(
                    f"'{param}' can't be None in a '{type(self).__name__}' "
                    "instance."
                )
            else:
                prelude = ""
                is_none = False

            if param == "min_atom_distance":
                comp = "greater than or equal to zero"
                valid = is_none or value >= 0
            else:
                if not is_none:
                    type_check(param, int)
                comp = "greater than zero"
                valid = is_none or value > 0
            msg = prelude + f"'{param}' must be {comp}, not {value}."
            if not valid:
                raise ValueError(msg)

        type_check("supports_slm_mask", bool)
        type_check("reusable_channels", bool)

        if not (0.0 < self.max_layout_filling <= 1.0):
            raise ValueError(
                "The maximum layout filling fraction must be "
                "greater than 0. and less than or equal to 1., "
                f"not {self.max_layout_filling}."
            )

        for ch_obj in self.channel_objects:
            type_check("All channels", GeneralChannel, value_override=ch_obj)
        

        if self.channel_ids is not None:
            if not (
                isinstance(self.channel_ids, (tuple, list))
                and all(isinstance(el, str) for el in self.channel_ids)
            ):
                raise TypeError(
                    "When defined, 'channel_ids' must be a tuple or a list "
                    "of strings."
                )
            if len(self.channel_ids) != len(set(self.channel_ids)):
                raise ValueError(
                    "When defined, 'channel_ids' can't have "
                    "repeated elements."
                )
            if len(self.channel_ids) != len(self.channel_objects):
                raise ValueError(
                    "When defined, the number of channel IDs must"
                    " match the number of channel objects."
                )
        else:
            # Make the channel IDs from the default IDs
            ids_counter: Counter = Counter()
            ids = []
            for ch_obj in self.channel_objects:
                id = ch_obj.default_id()
                ids_counter.update([id])
                if ids_counter[id] > 1:
                    # If there is more than one with the same ID
                    id += f"_{ids_counter[id]}"
                ids.append(id)
            object.__setattr__(self, "channel_ids", tuple(ids))

        if any(
            ch.basis == "XY" for ch in self.channel_objects
        ) and not isinstance(self.interaction_coeff_xy, float):
            raise TypeError(
                "When the device has a 'Microwave' channel, "
                "'interaction_coeff_xy' must be a 'float',"
                f" not '{type(self.interaction_coeff_xy)}'."
            )

        def to_tuple(obj: tuple | list) -> tuple:
            if isinstance(obj, (tuple, list)):
                obj = tuple(to_tuple(el) for el in obj)
            return obj

        # Turns mutable lists into immutable tuples
        for param in self._params():
            if "channel" in param:
                object.__setattr__(self, param, to_tuple(getattr(self, param)))

        self._initialize_couplings()
    def _validate_rydberg_states(self, ryd_stts:  tuple[tuple[int]], stt_lbls: tuple[str]|None) -> None:
        for state in ryd_stts:
            assert len(state)==4, f"Invalid Rydberg state: {state} - Rydberg state has to be defined with 4-length tuple, with parameters (n, l, j, ml)"

        if stt_lbls is None:
            self.state_labels = [ f"r{i}" for i,_ in enumerate(ryd_stts)]
        pass

    def _initialize_couplings(self) -> None:

        for i, (s1, lb1) in enumerate(zip(self.rydberg_states, self.state_labels)):
            for j, (s2, lb2) in enumerate(zip(self.rydberg_states[i:], self.state_labels[i:])):
                try:
                    c3 = getC3fromstates(self.alkali_atom, s1, s2, self.alkali_atom, s2, s1)
                except ValueError:
                    c3=0
                try:
                    c6 = getC6fromstates(self.alkali_atom, s1, self.alkali_atom, s2)
                except ValueError:
                    c6=0
                
                if c3!=0:
                    self.c3_dict[(lb1, lb2)] = c3
                    self.c3_dict[(lb2, lb1)] = c3
                if c6!=0:
                    self.c6_dict[(lb1, lb2)] = c6
                    self.c6_dict[(lb2, lb1)] = c6

                rvdw = getRvdw(self.alkali_atom, s1, self.alkali_atom, s2) # LeRoy radius in (\mu m)
                self.rvdw_dict[(lb1, lb2)] = rvdw
                self.rvdw_dict[(lb2, lb1)] = rvdw
        
    @property
    def _optional_parameters(self) -> tuple[str, ...]:
        return ()
    
    def _to_abstract_repr(self) -> dict[str, Any]:
        d = super()._to_abstract_repr()
        d["is_virtual"] = False
        return d
    def _to_dict(self) -> dict[str, Any]:
        return obj_to_dict(
            self, _build=False, _module="pulser.devices", _name=self.name
        )



