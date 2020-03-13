from setuptools import setup, find_packages

setup(
    name="ensembl_map",
    version="",
    author="mattdoug604",
    author_email="mattdoug604@gmail.com",
    packages=find_packages(),
    description="Map between gene, transcript, exon, CDS, and protein coordinates.",
    url="https://github.com/mattdoug604/ensembl_map.git",
    python_requires=">=3",
    install_requires=["pyensembl"],
)
