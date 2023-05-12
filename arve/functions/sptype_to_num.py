class sptype_to_num:

    def sptype_to_num(self, sptype:str) -> int:
        """Spectral type to number.

        :param sptype: spectral type
        :type sptype: str
        :return: spectral type represented as a number
        :rtype: int
        """

        return "OBAFGKM".index(sptype[0])*10 + int(sptype[1])