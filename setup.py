from setuptools import setup, find_packages

with open("requirements.txt", "r") as fin:
    requirements = list(map(str.strip, fin))

setup(
    name="multiseq",
    version="0.1.0",
    url="https://github.com/wflynny/multiseq-py",
    author="Bill Flynn",
    author_email="bill.flynn@jax.org",
    description="Python reimplementation of MULTI-seq deMULTIplex",
    packages=find_packages(),
    install_requires=requirements,
)
