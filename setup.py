from setuptools import setup, find_packages

setup(
    name="ensembl_map",
    version="",
    author="mattdoug604",
    author_email="mattdoug604@gmail.com",
    packages=find_packages(),
    description="Map between gene, transcript, exon, and protein coordinates.",
    url="",
    python_requires=">=3",
    install_requires=["pyensembl"],
)
