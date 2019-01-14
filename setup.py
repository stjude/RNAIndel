from setuptools import setup, find_packages

version = {}
with open("./RNAIndel/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="RNAIndel",
    version=version["__version__"],
    description="Somatic Indel Detector for Tumor RNA-Seq Data",
    url="https://github.com/adamdingliang/RNAIndel",
    author="Kohei Hagiwara, Liang Ding",
    author_email="kohei.hagiwara@stjude.org, liang.ding@stjude.org",
    license="Apache License 2.0",
    install_requires=[
        "pandas >= 0.23.0",
        "numpy >= 1.12.0",
        "scikit-learn == 0.18.1",
        "pysam == 0.13",
        "pyvcf == 0.6.8",
    ],
    python_requires=">=3.5.2",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    package_data={
        "RNAIndel": ["testdata/*"],
        "Bambino": [
            "bambino-1.0.jar",
            "mysql-connector-java-5.1.10.jar",
            "picard.jar",
            "third_party.jar",
            "testdata/*",
        ],
    },
    include_package_data=True,
    test_suite="tests",
    entry_points={
        "console_scripts": [
            "rna_indel=RNAIndel.rna_indel:main",
            "bambino=Bambino.bambino:main",
        ]
    },
)
