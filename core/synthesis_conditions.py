"""Optional synthesis metadata for match_peaks(). All fields are unused
by the current pipeline except `conditions_known`."""
from dataclasses import dataclass
from typing import Optional


@dataclass
class SynthesisConditions:
    activation_reagent:  Optional[str]   = None
    base:                Optional[str]   = None
    solvent:             Optional[str]   = None
    temperature_c:       Optional[float] = None
    method:              Optional[str]   = None
    coupling_time_min:   Optional[float] = None
    deprotection_cycles: Optional[int]   = None
    notes:               str             = ''

    def __post_init__(self):
        for attr in ('activation_reagent', 'base', 'solvent'):
            v = getattr(self, attr)
            if v:
                setattr(self, attr, v.strip())
        if self.method:
            self.method = self.method.lower().strip()

    @property
    def conditions_known(self):
        return any([self.activation_reagent, self.temperature_c, self.method])
