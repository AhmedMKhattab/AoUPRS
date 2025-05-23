from setuptools import setup, find_packages

# Read the full README for PyPI display
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='AoUPRS',
    version='0.2.0',
    description='AoUPRS is a Python module for calculating Polygenic Risk Scores (PRS) specific to the All of Us study',
    long_description=long_description,
    long_description_content_type='text/markdown',  # ← important
    author='Ahmed Khattab',
    packages=find_packages(),
    install_requires=[
        'hail',
        'gcsfs',
        'pandas',
    ],
    url='https://github.com/AhmedMKhattab/AoUPRS',  # optional but recommended
)