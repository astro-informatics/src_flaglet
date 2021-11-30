from pathlib import Path 
from skbuild import setup

cmake_args = ["-Dtests:BOOL=OFF"]

long_description = (
    Path(__file__).parent / "src" / "main" / "pyflaglet" / ".FLAGLET_Python_Documentation.rst"
).read_text()

setup(
    name="pyflaglet",
    version="1.0.0",
    author="Boris Leistedt & Jason McEwen & Matthew Price",
    install_requires=["numpy", "cython", "scipy"],
    extras_require={
        "dev": [
            "numpy"
            "setuptools",
            "wheel",
            "scikit-build",
            "cmake",
            "ninja",
            "cython",
            "conan",
            "black",
            "pytest-black"
            "pytest",
        ]
    },
    description="Fast wavelet transform on the ball",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/astro-informatics/src_flaglet",
    package_dir={"pyflaglet": "src/main/pyflaglet"},
    package_data={"pyflaglet": ["FLAGLET_Python_Documentation.rst"]},
    cmake_args=cmake_args,
    cmake_languages=("C",),
    license="GPL-3",
    packages=["pyflaglet"],
)
