from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='Metallicity_Stack_Commons',
    version='0.1.0',
    packages=find_packages('Metallicity_Stack_Commons'),
    url='https://github.com/astrochun/Metallicity_Stack_Commons',
    license='MIT License',
    author='Chun Ly',
    author_email='astro.chun@gmail.com',
    description='Set of common codes used in metallicity studies that use stacking techniques',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    install_requires = ['numpy', 'astropy', 'matplotlib', 'pylab', 'scipy']
)
