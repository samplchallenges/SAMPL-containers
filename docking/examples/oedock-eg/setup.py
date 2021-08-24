from setuptools import setup

setup(
    name='OE-Docking',
    version='0.1',
    py_modules=[
        "run_oedock", 
        "ligand", 
        "receptor",
        "docking",
	"oe_io"
    ],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        oedock=run_oedock:oedock
    ''',
)
