from setuptools import setup, find_packages

setup(
    name='ApChem',
    version='0.1.0',
    python_requires='>=3.7',
    package_dir={"":"."},
    packages=find_packages(".", include=["ApChem*",  "ApChem"]),
    install_requires=[
        'biopython',
        'freesasa==2.2.0a1',
        'gemmi',
        'numpy',
        'pandas',
        'pymp-pypi',
        'requests',
        'wget'
    ],
)