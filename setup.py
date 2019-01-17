from setuptools import setup, find_packages

version = {}
with open("./rnaindel/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="rnaindel",
    version=version["__version__"],
    description="Somatic indel discovery tool for tumor RNA-Seq data",
    url="https://github.com/stjude/RNAIndel",
    author="Kohei Hagiwara, Liang Ding",
    author_email="kohei.hagiwara@stjude.org, liang.ding@stjude.org",
    license="Apache License 2.0",
    install_requires=[
        "pandas >= 0.23.0",
        "numpy >= 1.12.0",
        "scikit-learn >= 0.18.1",
        "pysam >= 0.13",
    ],
    python_requires=">=3.5.2",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    entry_points={"console_scripts": ["rnaindel=rnaindel.rnaindel:main"]},
)
