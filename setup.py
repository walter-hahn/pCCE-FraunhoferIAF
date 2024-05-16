from setuptools import setup, find_packages

setup(
    name="pCCE",
    version="0.64",
    description="Calculation of the spin-coherence decay for quantum central-spin problems with focus on the Hahn echo problem.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/walter-hahn/pCCE-FraunhoferIAF",  # Replace with the actual URL
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v2.1",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
