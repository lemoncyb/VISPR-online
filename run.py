import yaml
import os
from vispr_screen import app
from vispr_screen.results import Screens
from vispr_screen.common import VisprError


def init_server(configs, host="127.0.0.1", port=5000):
    app.screens = Screens()
    print("Loading data.")
    print(configs)
    with open(configs) as f:
        config = yaml.load(f)
    try:
        app.screens.add(config, parentdir=os.path.dirname(configs))
    except KeyError as e:
        raise VisprError(
            "Syntax error in config file {}. Missing key {}.".format(configs,
                                                                         e))

yaml_config='/Users/cuiyb/workspace/PublicCRISPRScreenData/publiccrisprscreendata/24336569_LanderES_Science_2014/results/mle.vispr.yaml'
init_server(yaml_config)


app.run(debug=True)
