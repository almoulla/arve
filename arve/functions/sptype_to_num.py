class sptype_to_num:

    def sptype_to_num(
        self,
        sptype : str,
        ) -> int:
        """Spectral type to number.

        Parameters
        ----------
        sptype : str
            spectral type

        Returns
        -------
        int
            spectral type represented as a number
        """

        return "OBAFGKM".index(sptype[0])*10 + int(sptype[1])