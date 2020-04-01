from setuptools import setup, find_packages

name = "ensembl_map"
version = "1.0.2"

setup(
    name=name,
    version=version,
    author="mattdoug604",
    author_email="mattdoug604@gmail.com",
    packages=find_packages(),
    description="Convert between gene, transcript, exon, CDS, and protein coordinates using Ensembl.",
    url="https://github.com/mattdoug604/ensembl_map.git",
    python_requires=">=3",
    license="https://spdx.org/licenses/MIT.html",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=["pyensembl>=1.8.5"],
    tests_require=["mock", "nose"],
    test_suite="nose.collector",
)
