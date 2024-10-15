from setuptools import setup
import os


if os.path.isfile('requirements.txt'):
    with open('requirements.txt', 'r') as requirements_file:
        install_requires = requirements_file.read().splitlines()
for package in install_requires:
    if package.startswith('git'):
        pname = package.split('/')[-1].split('.')[0]
        install_requires[install_requires.index(package)] = pname + ' @ ' + package

setup(name='pyMALDIproc',
      version='0.4.6',
      url='https://github.com/gtluu/pyMALDIproc',
      license='Apache License',
      author='Gordon T. Luu',
      author_email='gtluu912@gmail.com',
      packages=['pymaldiproc', 'pymaldiviz', 'etc'],
      include_package_data=True,
      package_data={'': ['*.cfg']},
      description='pyMALDIproc: a Python package for working with MALDI mass spectrometry data',
      #entry_points={'console_scripts': ['pymaldiviz=pymaldiviz.launch_dashboard:main']},
      install_requires=install_requires,
      setup_requires=install_requires)
