from setuptools import setup

setup(
    name='AutoDock-singularity',
    version='0.1',
    py_modules=[
        'autodock',
        'main',
    ],
    install_requires=[
    ],
    entry_points='''
    [console_scripts]
    run-autodock=main:main_function
    '''
)

