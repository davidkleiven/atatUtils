
from setuptools import setup, find_packages

setup(
    name="atatutils",
    author=["David Kleiven"],
    author_email="davidkleiven4462gmail.com",
    long_description="Utility functions for atat",
    install_requires=['pymatgen', 'ase'],
    scripts=['scripts/traj2atat.py']
)