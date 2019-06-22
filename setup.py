
with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('find/version.py').read()) # loads __version__

setup(name='BGC_HMM_find',
      version=__version__,
      author='fernando garcia-guevara',
      description='Predict BGC presence in binned contigs',
      long_description=readme,
      license='GPL3+',
      keywords="",
      packages=find_packages(exclude='docs'),
      install_requires=('biopython >=1.64'),
      url='https://github.com/garciafertson/BGC_HMM_find',
      scripts=['bin/BGC_HMM_find.py', 'bin/BGC_HMM_train.py'],
)
                                    
