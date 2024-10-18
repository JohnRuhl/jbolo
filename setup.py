import os
from setuptools import find_packages, setup, Command


def readme():
    with open('README.md') as f:
        return f.read()

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf build dist *.pyc *.tgz *.egg-info')

setup(
    name="jbolo",
    package_dir={'jbolo': 'src'},
    packages=['jbolo'],    
    cmdclass={'clean':CleanCommand,},
    long_description = readme(),
    install_requires=[
        'numpy',
    ],
)
