import setuptools


def readme():
    with open("README.md", "r") as f:
        return f.read()


setuptools.setup(name="nbis-meta", version="0.0.1", author="John Sundh",
                 author_email="john.sundh@scilifelab.se",
                 description="A snakemake workflow for metagenomics",
                 long_description=readme(),
                 long_description_content_type="text/markdown",
                 license="MIT",
                 python_requires=">=3.7",
                 install_requires=[
                        "snakemake",
                    ],
                 url="https://github.com/nbisweden/nbis-meta",
                 packages=["workflow"],
                 package_data={"workflow": ["Snakefile", ".cfg/*", "envs/*",
                                            "notebooks/*", "report/*",
                                            "rules/*", "schemas/*",
                                            "scripts/*"]},
                 include_package_data=True,
                 entry_points={
                     'console_scripts': ['nmeta = workflow.__main__:main']},
                 classifiers=["Programming Language :: Python :: 3",
                              "License :: OSI Approved :: MIT License",
                              "Operating System :: OS Independent", ], )
