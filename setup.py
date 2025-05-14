from setuptools import setup

setup(
    name='reform',
    version='1.1.0',
    py_modules=['reform'],
    install_requires=[
        'biopython',
        'pgzip',
    ],
    entry_points={
        'console_scripts': [
            'reform = reform:main',
        ],
    },
    author="gencorefacility",
    description='A tool for editing genome sequences and annotation files',
    license='BSD-3-Clause',
)