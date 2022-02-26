
from modules.html_log.html_log import HtmlLog as log


log.config(filename="logs/my_test3.htm", mode='overwrite')


####################################################
class ClassOne:
    @staticmethod
    def class1_method1(s):
        log.write("Inside ClassOne.class1_method1()", indent=4)
        log.write(f"payload: {s}", indent=4)

####################################################


log.new_run()       # Log the start of a new run

log.write("Logging some stuff")

ClassOne.class1_method1("my data")

log.write("Logging more stuff")


log.separator("And, in summary:")   # Separator

log.write("Logging this in reverse font", style=log.reverse)
log.write("And now something blue, in italic", style=[log.italic, log.blue])