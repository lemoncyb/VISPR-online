from flask import Flask

app = Flask(__name__)
#app.jinja_env.trim_blocks = True
#app.jinja_env.lstrip_blocks = True

__version__ = '0.1.0'

from vispr import server