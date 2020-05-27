from setuptools import setup, find_namespace_packages


def readme():
    with open("README.md", "r") as f:
        return f.read()


setup(name="nbis-meta", version="0.0.0", author="John Sundh",
      author_email="john.sundh@scilifelab.se",
      description="A snakemake workflow for metagenomics",
      long_description=readme(), long_description_content_type="text/markdown",
      license="MIT", python_requires=">=3.7",
      install_requires=["snakemake", ],
      url="https://github.com/nbisweden/nbis-meta",
      packages=["nbis_meta"],
      package_dir={'nbis_meta': "workflow"},
      package_data={"nbis_meta": ["Snakefile", ".cfg/*", "envs/*", "notebooks/*",
                                 "report/*", "rules/*", "schemas/*",
                                 "scripts/*"]},
      entry_points={'console_scripts': ['nbis-meta = nbis_meta.__main__:main']},
      classifiers=["Programming Language :: Python :: 3",
                   "License :: OSI Approved :: MIT License",
                   "Operating System :: OS Independent", ], )
