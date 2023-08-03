from setuptools import setup, find_namespace_packages

setup(
    name='auto-routine-nanopore-qc',
    version='0.1.0-alpha',
    packages=find_namespace_packages(),
    entry_points={
        "console_scripts": [
            "auto-routine-nanopore-qc = auto_routine_nanopore_qc.__main__:main",
        ]
    },
    scripts=[],
    package_data={
    },
    install_requires=[
    ],
    description=' Automated routine quality control analysis of microbial sequence data from nanopore sequencers',
    url='https://github.com/BCCDC-PHL/auto-routine-nanopore-qc',
    author='Dan Fornika',
    author_email='dan.fornika@bccdc.ca',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
