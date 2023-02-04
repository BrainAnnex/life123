from src.modules.html_log.html_log import HtmlLog as log


import log_setup    # File that SETS UP some parameters of the logging


log.new_run()       # Log the start of a new run

log_setup.ClassOutside.method_x()

log.separator()
log.write("All is good!", style=log.bold)
