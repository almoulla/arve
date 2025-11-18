from setuptools import find_packages, setup

with open("README.md") as file:
    long_description = file.read()

setup(
    name="arve",
    version="0.6.0",
    author="Khaled Al Moulla",
    author_email="khaled.almoulla@gmail.com",
    description="Analyzing Radial Velocity Elements",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/almoulla/arve",
    license="MIT License",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=["astropy>=6.1.7"    ,
                      "astroquery>=0.4.10",
                      "lmfit>=1.3.3"      ,
                      "matplotlib>=3.10.1",
                      "numpy>=2.2.4"      ,
                      "pandas>=2.2.3"     ,
                      "scipy>=1.15.2"     ,
                      "tqdm>=4.67.1"      ],
    classifiers=["Development Status :: 1 - Planning"    ,
                 "Intended Audience :: Science/Research" ,
                 "License :: OSI Approved :: MIT License",
                 "Programming Language :: Python :: 3"   ],
    include_package_data=True,
    package_data={"arve": ["aux_data/masks/*.csv.zip",
                           "aux_data/spectra/*.csv.zip",
                           "aux_data/tellurics/*.csv.zip",
                           "aux_data/wavelengths/*.csv.zip"]}
)