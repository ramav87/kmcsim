from setuptools import setup, find_packages

with open("README.md", "r") as fr:
    long_description = fr.read()

setup(name='kmcsim',
      version='0.0.1',
      description='KMC simulation tool',
      long_description = long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/vlcekl/kmcsim',
      author='Lukas Vlcek',
      license='MIT',
      packages=find_packages(exclude='tests'),
      zip_safe=False)
