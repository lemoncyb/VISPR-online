import yaml
import os
from config import PUBLIC_DATA_DIR, PRIVATE_DATA_DIR, YAML_CONFIG
from vispr_screen import app
from vispr_screen.results import Screens
from vispr_screen.common import VisprError


def init_server(configs, condition, samples):
    app.screens = Screens()
    print("Loading data.")
    print(configs)
    with open(configs) as f:
        config = yaml.load(f)
    try:
        app.screens.add(config, condition, samples, parentdir=os.path.dirname(configs))
    except KeyError as e:
        raise VisprError(
            "Syntax error in config file {}. Missing key {}.".format(configs,
                                                                         e))


if __name__ == '__main__':

    dir_name = '24336569_LanderES_Science_2014'
    config_file = os.path.join(PUBLIC_DATA_DIR, dir_name, YAML_CONFIG)

    condition = 'HL60_final'  # used to select condition from gene_summary file
    initial_condition = 'HL60_initial,KBM7_initial'
    initial_conditions = [x.strip() for x in initial_condition.split(',')]

    if condition not in initial_conditions:
        samples = initial_conditions
        samples.append(condition) # used to select columns from count_normalized file
        init_server(config_file, condition, samples)
        app.run(debug=True)
    else:
        print condition + ' is initial condition!'