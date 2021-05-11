from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='Metallicity_Stack_Commons',
    version='1.4.7',
    packages=['Metallicity_Stack_Commons'],
    url='https://github.com/astrochun/Metallicity_Stack_Commons',
    license='MIT License',
    author='Chun Ly',
    author_email='astro.chun@gmail.com',
    description='Set of common codes used in metallicity studies that use stacking techniques',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=['numpy', 'astropy', 'matplotlib', 'scipy', 'requests',
                      'pytest', 'chun-codes']
)
