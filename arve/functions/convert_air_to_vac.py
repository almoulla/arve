import numpy as np


class convert_air_to_vac:
    def convert_air_to_vac(self, wave_air: np.ndarray) -> np.ndarray:
        """Convert wavelengths from air to vacuum.

        :param wave_air: wavelengths in air [Å]
        :return: wavelengths in vacuum [Å]
        """
        s = 1e4 / wave_air
        n = (
            1
            + 0.00008336624212083
            + 0.02408926869968 / (130.1065924522 - s**2)
            + 0.0001599740894897 / (38.92568793293 - s**2)
        )

        return wave_air * n
