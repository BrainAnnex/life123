from src.modules.html_log.html_log import HtmlLog as log


print("Prior to any setup, the (default) name for the log file is: ", log.LOG_FILENAME_BASE)

log.config(filename="logs/my_test2.htm", mode='overwrite')

print("After setup, the name for the log file is: ", log.LOG_FILENAME_BASE)

log.write("This is my log entry")

# The following line would give the error: "the config() method can only be invoked ONCE per run"
#log.config(filename="logs/my_test_new.htm", overwrite=True)
