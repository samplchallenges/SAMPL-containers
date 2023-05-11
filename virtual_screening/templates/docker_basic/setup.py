from setuptools import setup

setup(
    name='VirtualScreenContainer',
    version='0.1',
    py_modules=[
        'main',
    ],
    install_requires=[
    ],
    entry_points='''
        [console_scripts]
        run-virtual-screen=main:main_function
    '''
)
