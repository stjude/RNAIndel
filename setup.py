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
        "scikit-learn >= 0.18.1",
        "pysam >= 0.13.0", 
    ],
    python_requires=">=3.5.2",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    package_data={
        "rnaindel": [
            "bambino_lib/bambino-1.0.jar",
            "bambino_lib/mysql-connector-java-5.1.10.jar",
            "bambino_lib/picard.jar",
            "bambino_lib/third_party.jar",
        ]
    },
    entry_points={"console_scripts": ["rnaindel=rnaindel.rnaindel:main"]},
)
