from distutils.core import setup

setup(
    name='Contigo',
    version='0.1.0',
    author="Vince Forgetta",
    author_email="vincenzo.forgetta@mail.mcgill.ca",
    packages=['contigo'],
    scripts=['bin/contigo.py'],
    url='https://github.com/vforget/contigo',
    license='LICENCE.txt',
    description='Web-based assembly viewer.',
    long_description=open('README.txt').read(),
    install_requires=['PIL'],
)
