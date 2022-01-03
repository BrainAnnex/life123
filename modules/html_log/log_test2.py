from html_log import HtmlLog as log


# For testing purposes, we do the setup in this file
#log.config(filename="my_test.htm")
log.config(filename="my_test.htm", overwrite=True)
#log.config(filename="experiment_run.htm", multiple=True, max_files=3)

#log.config(filename="changed my mind.htm")


class ClassTwo:
    @staticmethod
    def class2_method1():
        log.write("Inside class2_method1()")
