from setuptools import find_packages, setup


setup(
    name="genefamilyflow",
    version="0.1.0",
    description="A configurable multi-species gene family analysis framework",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "click",
        "pyyaml",
        "biopython",
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "genefamilyflow=genefamilyflow.cli:main",
        ]
    },
)
