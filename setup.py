from setuptools import setup, find_packages

setup(
    name = "mds",
    version = "1.0",
    author = "Bingjia Yang",
    packages = find_packages(),
    install_requires=[
        'ase>=3.21',
        'numpy>=1.14.5',
    ]
)