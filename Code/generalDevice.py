from sys import version_info

from pulser.devices._device_datacls import Device, BaseDevice  
from dataclasses import dataclass, field, fields

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




@dataclass(frozen=True, repr=False)  # type: ignore[misc]
class GeneralDevice(BaseDevice):
    # Define C6 for arbitrary atoms
    # Compatible with the general framework
    name: str
    dimensions: DIMENSIONS
    alkali_atom: arc.AlkaliAtom = arc.Rubidium()
    rydberg_levels: tuple[tuple[int]] = None
    min_atom_distance: float
    max_atom_num: int | None
    max_radial_distance: int | None
    interaction_coeff_xy: float | None = None
    supports_slm_mask: bool = False
    max_layout_filling: float = 0.5
    reusable_channels: bool = field(default=False, init=False)
    channel_ids: tuple[str, ...] | None = None
    channel_objects: tuple[GeneralChannel, ...] = field(default_factory=tuple)
    pass


