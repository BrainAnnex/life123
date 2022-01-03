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


log.separator("Time for the plot:")

my_data = [
    {"id": 'd1', "index": 0, "value": 10, "region":'USA'},
    {"id": 'd2', "index": 1, "value": 11, "region":'Italy'},
    {"id": 'd3', "index": 2, "value": 12, "region":'Malta'},
    {"id": 'd4', "index": 3, "value": 6, "region":'Germany'}
]

handler = "D3_test1.js"

class_name = "svg1"

log.export_plot(class_name=class_name, js_file=handler, data=my_data)



exit(4)