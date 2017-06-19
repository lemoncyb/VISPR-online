from vispr_module.vispre.server import app

from vispr_module.vispre.cli import init_server


yaml_config='/Users/cuiyb/workspace/PublicCRISPRScreenData/publiccrisprscreendata/24336569_LanderES_Science_2014/results/mle.vispr.yaml'
init_server(yaml_config)

#app.run(debug=True)
