from setuptools import setup, find_packages

setup(
    name='Contigo',
    version='1.0.0',
    author="Vince Forgetta",
    author_email="vincenzo.forgetta@mail.mcgill.ca",
    packages=['contigo'],
    package_data={'contigo': ['contigo/static/*']},
    include_package_data=True,
    scripts=['bin/contigo'],
    url='https://github.com/vforget/contigo',
    license='LICENCE.txt',
    description='Web-based assembly viewer.',
    long_description=open('README').read(),
    requires=['PIL']
)
