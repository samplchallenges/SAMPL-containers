from setuptools import setup

setup(
    name='rdkit_logd',
    version='0.1',
    py_modules=['print_logd'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        print-LogD=print_logd:get_logd
    ''',
)
