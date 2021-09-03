from setuptools import setup

setup(
    name='containerName',
    version='0.1',
    py_modules=[
        'main_template',
    ],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        run-logd=main_template:logd_main
    '''
)
