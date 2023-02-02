"""
setup
"""

from setuptools import setup

setup(
    name="arve",
    version="0.0.2",
    description="Analyzing Radial Velocity Elements",
    url="https://github.com/almoulla/arve",
    author="Khaled Al Moulla",
    author_email="khaled.almoulla@gmail.com",
    license="MIT License",
    packages=["arve"],
    install_requires=["astroquery",
                      "lmfit"     ,
                      "matplotlib",
                      "numba"     ,
                      "numpy"     ,
                     ],

    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
)