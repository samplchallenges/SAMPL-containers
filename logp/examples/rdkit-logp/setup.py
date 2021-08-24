from setuptools import setup

setup(
    name='rdkit_logP',
    version='0.1',
    py_modules=['print_logp'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        print-LogP=print_logp:get_logp
    ''',
)
