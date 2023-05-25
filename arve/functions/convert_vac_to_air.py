import numpy as np


class convert_vac_to_air:
    def convert_vac_to_air(self, wave_vac: np.ndarray) -> np.ndarray:
        """Convert wavelengths from vacuum to air.

        :param wave_vac: wavelengths in vacuum [Å]
        :return: wavelengths in air [Å]
        """
        s = 1e4 / wave_vac
        n = (
            1
            + 0.0000834254
            + 0.02406147 / (130 - s**2)
            + 0.00015998 / (38.9 - s**2)
        )

        return wave_vac / n
