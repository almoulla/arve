"""
add vpsd component
"""


def add_vpsd_component(star, label, name, type, coef, vary):

    # add label if not already existing
    if label not in star.vpsd_components:
        star.vpsd_components[label] = {}

    # save VPSD component
    star.vpsd_components[label][name] = {"type": type, "coef": coef, "vary": vary}