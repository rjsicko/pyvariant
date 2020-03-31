from setuptools import setup, find_packages

setup(
    name="ensembl_map",
    version="1.0.0",
    author="mattdoug604",
    author_email="mattdoug604@gmail.com",
    packages=find_packages(),
    description="Convert between gene, transcript, exon, CDS, and protein coordinates using Ensembl.",
    url="https://github.com/mattdoug604/ensembl_map.git",
    python_requires=">=3.3",
    install_requires=["pyensembl>=1.8.5"],
    tests_require=["mock", "nose"],
    test_suite="nose.collector",
)
