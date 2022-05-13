from setuptools import find_packages, setup

setup(
    name="ensembl_map",
    version="2.0.0",
    author="Matthew Douglas",
    author_email="mattdoug604@gmail.com",
    description="Map sequence variants between chromosome, gene, exon, protein, and transcript representations.",
    url="https://github.com/mattdoug604/ensembl_map.git",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "BioPython>=1.73",
        "logzero",
        "pyensembl",
        "pyfaidx",
    ],
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
)
