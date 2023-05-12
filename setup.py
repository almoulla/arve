from setuptools import find_packages, setup

with open("README.md", "r") as file:
    long_description = file.read()

setup(
    name="arve",
    version="0.1.3",
    author="Khaled Al Moulla",
    author_email="khaled.almoulla@gmail.com",
    description="Analyzing Radial Velocity Elements",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/almoulla/arve",
    license="MIT License",
    packages=find_packages(),
    install_requires=["astroquery",
                      "lmfit"     ,
                      "matplotlib",
                      "numba"     ,
                      "numpy"     ,
                      "scipy"     ,
                      "tqdm"      ,
                     ],
    classifiers=["Development Status :: 1 - Planning"    ,
                 "Intended Audience :: Science/Research" ,
                 "License :: OSI Approved :: MIT License",
                 "Programming Language :: Python :: 3"   ,
                 ],
    include_package_data=True,
    package_data={"arve": ["aux_data/masks/*.csv"]}
)