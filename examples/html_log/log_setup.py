from life123.html_log import HtmlLog as log

# This file sets up the log configuration, and also provides a class


log.config(filename="logs/test_with_ext_setup.htm", mode='overwrite')

#log.config(filename="experiment_run.htm", mode='multiple', max_files=3)
#log.config(filename="my_test.htm")
#log.config(filename="changed my mind.htm")


class ClassOutside:
    @staticmethod
    def method_x():
        log.write("Inside ClassOutside.method_x()", indent=4)

