from setuptools import find_packages, setup

setup(
    name="ensembl_map",
    version="2.0.0",
    author="Matthew Douglas",
    author_email="mattdoug604@gmail.com",
    description="Convert between cDNA, DNA, exon, protein, and RNA coordinates.",
    url="https://github.com/mattdoug604/ensembl_map.git",
    packages=find_packages(),
    include_package_data=True,
    python_requires=">=3.8",
    install_requires=["BioPython>=1.73", "pyfaidx"],
    extras_require={
        "dev": [
            "black",
            "flake8",
            "isort",
            "mypy",
            "pytest",
            "pytest-cov",
            "pytest-mock",
            "twine",
            "wheel",
        ]
    },
    entry_points={"console_scripts": ["ensembl-map = ensembl_map.cli:main"]},
)
