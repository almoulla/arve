from .doppler_shift     import doppler_shift
from .gls_periodogram   import gls_periodogram
from .inverted_gaussian import inverted_gaussian
from .sptype_to_num     import sptype_to_num

class _Functions_classes(doppler_shift,
                         gls_periodogram,
                         inverted_gaussian,
                         sptype_to_num):
    pass