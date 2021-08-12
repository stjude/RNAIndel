from setuptools import setup, find_packages

version = {}
with open("./rnaindel/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="rnaindel",
    version=version["__version__"],
    description="Somatic indel discovery tool for tumor RNA-Seq data",
    url="https://github.com/stjude/RNAIndel",
    author="Kohei Hagiwara",
    author_email="kohei.hagiwara@stjude.org",
    license="Apache License 2.0",
    install_requires=[
        "pandas >= 0.23.0",
        "scikit-learn >= 0.22.0",
        "indelpost >= 0.0.4", 
    ],
    python_requires=">=3.6",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    package_data={
        "rnaindel": [
            "defaultcaller/bambino-1.0.jar",
            "defaultcaller/mysql-connector-java-5.1.10.jar",
            "defaultcaller/picard.jar",
            "defaultcaller/third_party.jar",
        ]
    },
    entry_points={"console_scripts": ["rnaindel=rnaindel.rnaindel:main"]},
)
