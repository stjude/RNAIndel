from setuptools import setup, find_packages

version = {}
with open("./rna_indel/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="RNAIndel",
    version=version["__version__"],
    description="Somatic indel discovery tool for tumor RNA-Seq data",
    url="https://github.com/stjude/RNAIndel",
    author="Kohei Hagiwara, Liang Ding",
    author_email="kohei.hagiwara@stjude.org, liang.ding@stjude.org",
    license="Apache License 2.0",
    install_requires=[
        "pandas >= 0.23.0",
        "numpy >= 1.12.0",
        "scikit-learn == 0.18.1",
        "pysam >= 0.13",
    ],
    python_requires=">=3.5.2",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    package_data={
        "rna_indel": ["bambino_lib/bambino-1.0.jar", "bambino_lib/mysql-connector-java-5.1.10.jar", "bambino_lib/picard.jar", "bambino_lib/third_party.jar"]
    },
    include_package_data=True,
    test_suite="tests",
    # I want have only one entry ponint
    entry_points={
        "console_scripts": [
            "rna_indel=rna_indel.rna_indel:main",
        ]
    },
)
