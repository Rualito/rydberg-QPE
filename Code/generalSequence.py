import pulser.sampler as sampler
from pulser import Pulse, Sequence
from pulser.sequence._call import _Call

import numpy as np

from typing import Any, Optional, Tuple, Union, cast, overload
from pulser.register.base_register import BaseRegister, QubitId
from collections.abc import Iterable, Mapping

from generalDevice import GeneralDevice
from pulser.register.mappable_reg import MappableRegister


class GeneralSequence(Sequence):
    pass
    # Methods to overload:
    # available_channels()
    # magnetic_field()
    # set_magnetic_field()
    # declare_channel()
    # measure()
    def __init__(self, register: Union[BaseRegister, MappableRegister], device: GeneralDevice):
        super().__init__(register, device)
    def declare_channel(self, name: str, channel_id: str, initial_target: Optional[Union[QubitId, Iterable[QubitId]]] = None) -> None:
        self._in_multi = True
        self.set_magnetic_field_multi()
        
        return super().declare_channel(name, channel_id, initial_target)
    

    def set_magnetic_field_multi(self, bx: float = 0, by: float = 0, bz: float = 30) -> None:
        """Sets the magnetic field acting on the entire array. Multi-Rydberg mode.

        The magnetic field vector is defined on the reference frame of the
        atoms in the Register (with the z-axis coming outside of the plane).
        Can only be defined before there are pulses added to the sequence.

        Note:
            The magnetic field only work in the "XY Mode". If not already
            defined through the declaration of a Microwave channel, calling
            this function will enable the "XY Mode".

        Args:
            bx: The magnetic field in the x direction (in Gauss).
            by: The magnetic field in the y direction (in Gauss).
            bz: The magnetic field in the z direction (in Gauss).
        """
        if not self._in_multi:
            if self._schedule:
                raise ValueError(
                    "The magnetic field can only be set in 'XY Mode'."
                )
            # No channels declared yet
            # self._in_xy = True
            self._in_multi = True
        elif not self._empty_sequence:
            # Not all channels are empty
            raise ValueError(
                "The magnetic field can only be set on an empty sequence."
            )

        mag_vector = (bx, by, bz)
        if np.linalg.norm(mag_vector) == 0.0:
            raise ValueError(
                "The magnetic field must have a magnitude greater than 0."
            )
        self._mag_field = mag_vector

        # No parametrization -> Always stored as a regular call
        self._calls.append(_Call("set_magnetic_field", mag_vector, {}))

