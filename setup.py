from distutils.core import setup

setup(
    name='chimera_autoguider',
    version='0.0.1',
    packages=['chimera_autoguider', 'chimera_autoguider.instruments', 'chimera_autoguider.controllers'],
    scripts=['scripts/chimera-autoguide'],
    url='http://github.com/astroufsc/chimera_autoguider',
    license='GPL v2',
    author='Tiago Ribeiro',
    author_email='tribeiro@ufs.br',
    description='Chimera controller plugin for autoguiding telescopes'
)
