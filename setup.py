from setuptools import setup, find_namespace_packages


def readme():
    with open("README.md", "r") as f:
        return f.read()


setup(name="nbismeta", version="0.0.0", author="John Sundh",
      author_email="john.sundh@scilifelab.se",
      description="A snakemake workflow for metagenomics",
      long_description=readme(), long_description_content_type="text/markdown",
      license="MIT", python_requires=">=3.7",
      install_requires=["snakemake", ],
      url="https://github.com/nbisweden/nbis-meta",
      packages=["nbismeta"],
      package_data={"nbismeta": ["Snakefile", ".cfg/*", "envs/*", "notebooks/*",
                                 "report/*", "rules/*", "schemas/*",
                                 "scripts/*"]},
      entry_points={'console_scripts': ['nbismeta = nbismeta.__main__:main']},
      classifiers=["Programming Language :: Python :: 3",
                   "License :: OSI Approved :: MIT License",
                   "Operating System :: OS Independent", ], )
