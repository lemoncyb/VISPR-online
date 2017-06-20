# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

__version__ = '0.1.0'


from flask import Flask

app = Flask(__name__)
#app.jinja_env.trim_blocks = True
#app.jinja_env.lstrip_blocks = True

from vispr_screen import server

