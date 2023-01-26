"""
add vpsd component
"""


def add_vpsd_component(star, name, type, coef, vary):

    # save VPSD component
    star.vpsd_components[name] = {"type": type, "coef": coef, "vary": vary}