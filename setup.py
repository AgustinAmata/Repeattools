from setuptools import setup
from os.path import join, dirname
setup(
    name='Repeattools',
    version='1.0',
    packages=["Repeattools"],
    long_description=open(join(dirname(__file__), 'README.md')).read(),
)