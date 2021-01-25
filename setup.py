import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="grubin",
    version="0.0.1",
    author="Aaron Johnson",
    author_email="johnsoad@uwm.edu",
    description="Gelman-Rubin Rhat statistic for use with PTMCMCSampler output.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AaronDJohnson/grubin",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)