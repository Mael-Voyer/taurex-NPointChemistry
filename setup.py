from setuptools import setup, find_packages

packages = find_packages(exclude=("tests", "doc"))

provides = ["taurex_npoint_chemistry"]

install_requires = [
    "taurex",
]

entry_points = {
    "taurex.plugins": [
        "npoint_chemistry = taurex_npoint_chemistry",
    ]
}

setup(
    name="taurex-npoint-chemistry",
    author="Mael Voyer",
    version="1.0",
    author_email="mael.voyer@cea.fr",
    license="BSD-3-Clause",
    description="N-point free chemistry profile for TauREx",
    packages=packages,
    entry_points=entry_points,
    provides=provides,
    install_requires=install_requires,
)