from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='Metallicity_Stack_Commons',
    version='1.2.0',
    packages=['Metallicity_Stack_Commons'],
    url='https://github.com/astrochun/Metallicity_Stack_Commons',
    license='MIT License',
    author='Chun Ly',
    author_email='astro.chun@gmail.com',
    description='Set of common codes used in metallicity studies that use stacking techniques',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'numpy>=1.18.4',
        'astropy>=3.2.2',
        'matplotlib>=3.1.1',
        'scipy>=1.3.1',
        'requests>=2.22.0',
        'pytest'
    ],
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ]
)
