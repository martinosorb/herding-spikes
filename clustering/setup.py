import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
#def read(fname):
#    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "herding-spikes",
    version = "0.0.4",
    author = "Martino Sorbaro",
    author_email = "martino.sorbaro@ed.ac.uk",
    description = "Software for high density multi-electrode array recordings -- clustering",
    license = "GPL",
    keywords = "spikes neuro sorting",
    url = "https://github.com/martinosorb/herding-spikes",
    py_modules=['herdingspikes','ms_new_withmp'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: GPL License",
    ],
)
