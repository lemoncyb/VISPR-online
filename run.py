from invoke_vispr.vispre.server import app

from invoke_vispr.vispre.cli import init_server


config='/Users/cuiyb/workspace/PublicCRISPRScreenData/publiccrisprscreendata/24336569_LanderES_Science_2014/results/mle.vispr.yaml'
init_server(config)

#app.run(debug=True)
