import setuptools
def readme():
    with open("README.md", "r") as f:
        return f.read()

setuptools.setup(
    name="nbis-meta",
    version="0.0.1",
    author="John Sundh",
    author_email="john.sundh@scilifelab.se",
    description="A snakemake workflow for metagenomics",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/nbisweden/nbis-meta",
    package_dir={"": "."},
    packages=["workflow"],
    package_data={"workflow": ["workflow/*", "config/*"]},
    entry_points={'console_scripts': ['workflow = workflow.__main__:main']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
