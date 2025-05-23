from setuptools import setup, find_packages

setup(
    name='AoUPRS',
    version='0.1.5',
    description='AoUPRS is a Python module for calculating Polygenic Risk Scores (PRS) specific to the All of Us study',
    author='Ahmed Khattab',
    packages=find_packages(),
    install_requires=[
        'hail',
        'gcsfs',
        'pandas',
    ],
    long_description='AoUPRS is a Python module for calculating Polygenic Risk Scores (PRS) specific to the All of Us study. It integrates with Hail for genetic data processing.'
)
