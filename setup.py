from setuptools import setup
from os.path import join, dirname
setup(
    name='Repeattools',
    version='1.0',
    author="Agustin Amata",
    author_email="agus.amata2002@gmail.com",
    url="https://github.com/AgustinAmata/Repeattools.git",
    packages=["Repeattools"],
    long_description=open(join(dirname(__file__), 'README.md')).read(),
    install_requires=["ete3==3.1.3", "matplotlib==3.7.3",
                      "numpy==1.25.0", "pandas==2.0.3",
                      "scipy==1.11.1", "seaborn==0.12.2",
                      "scikit-learn==1.3.0"]
)
