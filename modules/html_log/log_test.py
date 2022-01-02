from html_log import HtmlLog as log
print(log.LOG_FILENAME_BASE)

import log_test2

#log.config(filename="another catchy name.htm")

class ClassOne:
    @staticmethod
    def class1_method1(s):
        log.write("Inside class1_method1()")
        log.write(f"payload: {s}")

####################################################

log.new_run()

print(log.LOG_FILENAME_BASE)

log_test2.ClassTwo.class2_method1()

log.write("Alles ist gut!")

ClassOne.class1_method1("my 5th round of data")

log.separator("And, in summary:")

#log.write("Wunderbar!", style=[log.italic, log.blue])
log.write("Wunderbar!", style=log.reverse)

#log.config(filename="another catchy name.htm")



exit(4)