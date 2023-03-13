
from pulser.channels.base_channel import Channel
import arc
from dataclasses import dataclass, field, fields

from sys import version_info
if version_info[:2] >= (3, 8):  # pragma: no cover
    from typing import Literal
else:  # pragma: no cover
    try:
        from typing_extensions import Literal  # type: ignore
    except ImportError:
        raise ImportError(
            "Using pulser with Python version 3.7 requires the"
            " `typing_extensions` module. Install it by running"
            " `pip install typing-extensions`."
        )

import numpy as np

default_local_channel_params ={
    'max_abs_detuning':2 * np.pi * 20,
    'max_amp':2 * np.pi * 10,
    'min_retarget_interval':220,
    'fixed_retarget_t':0,
    'max_targets':1,
    'clock_period':4,
    'min_duration':16,
    'max_duration':2**26,
}

default_global_channel_params ={
    'max_abs_detuning':2 * np.pi * 20,
    'max_amp':2 * np.pi * 10,
    'clock_period':4,
    'min_duration':16,
    'max_duration':2**26,
}

@dataclass(init=True, repr=False, frozen=True)  # type: ignore[misc]
class GeneralChannel(Channel):
    """
    <Insert description of class> 
    """
    alkali_atom: arc.AlkaliAtom=arc.Rubidium()
    state1_lbl: str = 'r1'
    state2_lbl: str = 'r2'
    basis: str = ''  # arbitrary name for the basis, only has to identify transition
    state_dict: dict = field(default_factory=dict)
    
    
    @property
    def _internal_param_valid_options(self) -> dict[str, tuple[str, ...]]:
        """Internal parameters and their valid options."""
        return dict(
            name=("Rydberg", "Raman", "Microwave", 'GeneralChannel'),
            basis=("ground-rydberg", "digital", "XY", self.basis),
            addressing=("Local", "Global"),
        )

    def __post_init__(self) -> None:
        """Validates the channel's parameters."""
        super().__post_init__()

        if ( self.state1_lbl not in ('g','h') ) and ( self.state2_lbl not in ('g', 'h')): 
            GeneralChannel.validate_transition(self.alkali_atom, self.state_dict[self.state1_lbl], self.state_dict[self.state2_lbl])

    
    def validate_transition(atom:arc.AlkaliAtom, state1:tuple[int], state2:tuple[int],):
        """
        Checks if the transition state1 <-> state2 is valid through field excitation.
        """ 
        # w: 1/e^2 beam radius
        # P: laser power
        # atom.getRabiFrequency(n1,l1,j1,mj1,n2,l2,j2 ,mj2,P,w)

        assert np.abs(atom.getRabiFrequency(*state1,*state2,1,1))>1e-8, "Invalid transition: Rabi frequency is 0. Check selection rules (orbital momentum l or total momentum j)."

    # @property
    def basis(self) -> str:
        """The addressed basis name."""
        return str((self.state1_lbl, self.state2_lbl))