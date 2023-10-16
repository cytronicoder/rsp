from setuptools import setup, find_packages

with open("requirements.txt") as f:
    required = f.read().splitlines()

setup(
    name="rsp",
    version="0.1",
    packages=find_packages(),
    install_requires=required,
    author="Zeyu Yao",
    author_email="novodoodle@gmail.com",
    description="A novel method offering visual representation based on cell coverage and distribution within a high-dimensional spatial embedding space.",
    license="MIT",
    keywords="visualization cell coverage distribution t-sne rsp",
)
